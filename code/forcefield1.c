// file forcefield1.c
// invoke methods in following order:
//   FF_new, FF_set_<parm>, FF_build
// then invoke following in any order:
//   FF_get_<parm>, FF_energy, FF_rebuild
// finally invoke FF_delete
#include <assert.h>  //:::
#include <stdio.h>  //: might remove later
#include <stdlib.h>
#include <stdbool.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "forcefield.h"



static void neighborlist(FF *ff, int N, Vector *position,double margin);

// particles are not assumed to be in periodic cell.
// helper functions
static double partcl2partcl(FF *ff, int N, Vector *force, Vector *position,
                            double *charge);
static void anterpolate(FF *ff, Triple gd, double *q, int N, double *charge,
                        Vector *r);
static void restrict_(FF *ff, Triple gd, double *ql, double *qlm1, int lm1);
static void prolongate(FF *ff, Triple gd, double *el, double *elp1);
void grid2grid(FF* restrict ff,const int l,const Triple gd,
               double* restrict el, const double* restrict ql,
               const Triple sd, const double* restrict kh);
static void FFTg2g(FF *ff, Triple gd, double *el, double *ql, double *kh);
static void DFTg2g(Triple gd, double *el, double *ql, double *kh);
static void interpolate(FF *ff, int N, Vector *E, Vector *r, Triple gd,
                          double *el);
double msm4g_tictocmanager(int push) {
  double elapsed_seconds = 0.0;
  static clock_t stack_data[100] ;
  static int stack_lastindex = 0 ;
  if (push) {
    stack_data[stack_lastindex] = clock();
    stack_lastindex = stack_lastindex + 1;
  } else {
    clock_t now = clock();
    stack_lastindex = stack_lastindex - 1;
    clock_t previous = stack_data[stack_lastindex];
    elapsed_seconds = (double)(now-previous)/CLOCKS_PER_SEC;
  }
  return elapsed_seconds;
}
void msm4g_tic() {
  msm4g_tictocmanager(1);
}

double msm4g_toc() {
  return msm4g_tictocmanager(0);
}

void padding(FF *ff,double* ql, Triple gd, int l){
  msm4g_tic();
  int nu = ff->orderAcc;
  int dmax = 2*ff->nLim + 1;
  Triple sd, of;
  sd.x = max(min(dmax,gd.x),nu);
  sd.y = max(min(dmax,gd.y),nu);
  sd.z = max(min(dmax,gd.z),nu);
  of.x = max(sd.x/2, nu/2);
  of.y = max(sd.y/2, nu/2);
  of.z = max(sd.z/2, nu/2);
  
  int gdznew = gd.z + sd.z;
  int gdynew = gd.y + sd.y;
  int gdxnew = gd.x + sd.x;
  
  for (int i = 0 ; i < gdxnew ;i++) {
    int qx = i;
    while (qx < of.x)            qx += gd.x;
    while (qx > gd.x - 1 + of.x) qx -= gd.x;
    for (int j = 0 ; j < gdynew ;j++) {
      int qy = j;
      while (qy < of.y)            qy += gd.y;
      while (qy > gd.y - 1 + of.y) qy -= gd.y;
      for (int k = 0 ; k < gdznew ;k++) {
        int qz = k;
        while (qz < of.z)            qz += gd.z;
        while (qz > gd.z - 1 + of.z) qz -= gd.z;
        
        if (i < of.x || i > of.x + gd.x - 1 ||
            j < of.y || j > of.y + gd.y - 1 ||
            k < of.z || k > of.z + gd.z - 1 ) {
          int inside  = qx * gdynew * gdznew + qy * gdznew + qz;
          int outside =  i * gdynew * gdznew +  j * gdznew + k;
          
          ql[outside] = ql[inside];
        }
      }
    }
  }
  ff->time_padding[l] = msm4g_toc();
}

double FF_energy(FF *ff, double (*force)[3], double (*position)[3], double *weight) {
  // N may change in case of grand canonical simulations
  double *tau = ff->tau;
  double *charge = ff->q;
  int nu = ff->orderAcc;
  int N = ff->N;
  Vector *r = (Vector *)position;
  Vector *F = (Vector *)force;
  double *wt = (double *)malloc((ff->maxLevel + 1)*sizeof(double));
  for (int l = 0; l <= ff->maxLevel; l++)
    wt[l] = weight == NULL ? 1. : weight[l];
  // compute forces and energy
  // calculate the constant part of energy
  double qsum2 = 0, q2sum = 0.;
  for (int i = 0; i < N; i++)
    qsum2 += charge[i], q2sum += charge[i]*charge[i];
  qsum2 *= qsum2;
  //printf("csr: %25.16f\n",0.5*qsum2*ff->coeff2);
  //printf("ushort_self:: %25.16f\n",0.5*q2sum*ff->coeff1);
  double energy = 0.5*q2sum*ff->coeff1 - 0.5*qsum2*ff->coeff2;
  
  
  double s = 0 - 1.;
  // gamma(0) = tau(0 - 1)
  double gam = tau[nu-1];
  for (int k = nu - 2; k >= 0; k--)
    gam = tau[k] + s*gam;
  double ulong_self = - 0.5*q2sum * gam/ff->aCut[0];
  //printf("ulong_self: %25.16f\n",ulong_self);
  energy += ulong_self ;
  // particle-to-particle
  energy *= wt[0];
  msm4g_tic();
  double ushort_real = partcl2partcl(ff, N, F, r, charge) ;
  ff->time_partcl2partcl = msm4g_toc();
  //printf("ushort_real : %25.16f\n",ushort_real);
  energy += wt[0]*ushort_real;
  for (int i = 0; i < N; i++)
    F[i].x *= wt[0], F[i].y *= wt[0], F[i].z *= wt[0];
  double **q = (double **)malloc((ff->maxLevel + 1)*sizeof(double *));
  Triple gd = *(Triple *)ff->topGridDim;
  // no padding at the top-level
  q[ff->maxLevel] = (double *)calloc((gd.x+nu-1)*(gd.y+nu-1)*(gd.z+nu-1), sizeof(double));
  for (int l = ff->maxLevel-1; l >= 1; l--){
    gd.x *= 2; gd.y *= 2; gd.z *= 2;
    int dmax = 2*ff->nLim + 1;
    Triple sd = {max(min(dmax,gd.x),nu-1),
                 max(min(dmax,gd.y),nu-1),
                 max(min(dmax,gd.z),nu-1)};
    
    q[l] = (double *)calloc((gd.x+sd.x)*(gd.y+sd.y)*(gd.z+sd.z), sizeof(double));
    }
  // calculate q[1]
  msm4g_tic();
  anterpolate(ff, gd, q[1], N, charge, r);
  ff->time_anterpolation = msm4g_toc();
  for (int l = 2; l <= ff->maxLevel; l++){
    padding(ff, q[l-1], gd, l-1);
    gd.x /= 2; gd.y /= 2; gd.z /= 2;
    // calculate q[l]
    msm4g_tic();
    restrict_(ff, gd, q[l], q[l-1], l-1);
    ff->time_restriction[l]= msm4g_toc();
  }
  double *el = (double *)calloc(gd.x*gd.y*gd.z, sizeof(double));
  // set e^L = calK^L q^L
  Triple sd = gd;
  int l = ff->maxLevel;
  for (int m = 0; m < gd.x*gd.y*gd.z; m ++) q[l][m] *= wt[l];
  msm4g_tic();
    if (ff->FFT) FFTg2g(ff, gd, el, q[l], ff->khat[l]);
    else DFTg2g(gd, el, q[l], ff->khat[l]);
  ff->time_grid2grid[l] = msm4g_toc();
  for (int l = ff->maxLevel-1; l >= 1; l--){
    free(q[l + 1]);
    double *elp1 = el;
    gd.x  *= 2; gd.y *= 2;  gd.z *= 2;
    el = (double *)calloc(gd.x*gd.y*gd.z, sizeof(double));
    // add prolongated e^{l+1} to e^l
    msm4g_tic();
    prolongate(ff, gd, el, elp1);
    ff->time_prolongation[l] = msm4g_toc();
    free(elp1);
    for (int m = 0; m < gd.x*gd.y*gd.z; m ++) q[l][m] *= wt[l];
    // add calK^l q^l to e^1
    int dmax = 2*ff->nLim + 1;
    Triple sd = {dmax < gd.x ? dmax : gd.x,
                 dmax < gd.y ? dmax : gd.y,
                 dmax < gd.z ? dmax : gd.z};
    
    grid2grid(ff,l,gd, el, q[l], sd, ff->khat[l]);
    }
  free(wt);
  // add in grid-level energy
  int dmax = 2 * ff->nLim + 1;
  int padx = 0, pady = 0, padz = 0;
  if (ff->maxLevel > 1) {
    padx = max(min(dmax,gd.x),nu);
    pady = max(min(dmax,gd.y),nu);
    padz = max(min(dmax,gd.z),nu);
  }
  int ofx = padx / 2 ;
  int ofy = pady / 2 ;
  int ofz = padz / 2 ;

  double grid_en = 0.;
  for (int i = 0 ; i < gd.x; i++) {
    for (int j = 0 ; j < gd.y; j++) {
      for (int k = 0; k < gd.z; k++) {
        int elidx = (i * gd.y + j) * gd.z + k;
        int q1idx = ((i + ofx) * (gd.y + pady) + (j+ofy)) * (gd.z + padz) + (k + ofz);
        grid_en += q[1][q1idx]*el[elidx];
      }
    }
  }
  double ulong = 0.5*grid_en ;
  //printf("ulong_direct + fourier : %25.16e\n",ulong);
  
  energy += 0.5*grid_en;
  free(q[1]);
  free(q);
  // initialize grid-level electric field
  Vector *E = (Vector *)malloc(N*sizeof(Vector));
  msm4g_tic();
  interpolate(ff, N, E, r, gd, el);
  ff->time_interpolation = msm4g_toc();
  free(el);
  for (int i = 0; i < N; i++){
    F[i].x += charge[i]*E[i].x;
    F[i].y += charge[i]*E[i].y;
    F[i].z += charge[i]*E[i].z;}
  free(E);
  ff->errEst = FF_get_errEst(ff);
  return energy;
}




static double partcl2partcl(FF *ff, int N, Vector *force, Vector *position,
                            double *charge){
  double energy = 0.0;
  int nu = ff->orderAcc;
  double a_0 = ff->aCut[0];
  double a_02 = a_0 * a_0 ;
  Matrix A = *(Matrix *)ff->A;
  Matrix Ai = *(Matrix *)ff->Ai;
  for (int i = 0; i < N; i++) force[i].z = force[i].y = force[i].x = 0.;
  Vector as = {sqrt(Ai.xx*Ai.xx + Ai.yx*Ai.yx + Ai.zx*Ai.zx),
               sqrt(Ai.xy*Ai.xy + Ai.yy*Ai.yy + Ai.zy*Ai.zy),
               sqrt(Ai.xz*Ai.xz + Ai.yz*Ai.yz + Ai.zz*Ai.zz)};
  double pxlim = ceil(a_0*as.x - 0.5);
  double pylim = ceil(a_0*as.y - 0.5);
  double pzlim = ceil(a_0*as.z - 0.5);
 
  
  int *nlist = ff->nlist;
  int ptr = 0;
  int i = nlist[ptr++];
  int int_list = 0;
  while (i != -1) {
    Vector ri = position[i];
    int j = nlist[ptr++];
    while ( j != -1 ) {
      Vector rj = position[j];
      Vector r = {rj.x - ri.x, rj.y - ri.y, rj.z - ri.z};
      // convert to nearest image
      //Vector s = prod(Ai, r);
      Vector s = {Ai.xx*r.x + Ai.xy*r.y + Ai.xz*r.z,
                  Ai.yx*r.x + Ai.yy*r.y + Ai.yz*r.z,
                  Ai.zx*r.x + Ai.zy*r.y + Ai.zz*r.z};
 
      Vector p = {floor(s.x + 0.5), floor(s.y + 0.5), floor(s.z + 0.5)};
      //Vector Ap = prod(A, p);
      Vector Ap ={A.xx*p.x + A.xy*p.y + A.xz*p.z,
                  A.yx*p.x + A.yy*p.y + A.yz*p.z,
                  A.zx*p.x + A.zy*p.y + A.zz*p.z};
      r.x -= Ap.x; r.y -= Ap.y; r.z -= Ap.z;
      Vector r0 = r;
      double sum = 0.;
      Vector F = {0., 0., 0.};
      int pcount = 0;
      for (p.x = - pxlim; p.x <= pxlim; p.x++)
        for (p.y = - pylim; p.y <= pylim; p.y++)
          for (p.z = - pzlim; p.z <= pzlim; p.z++){
            //Ap = prod(A, p);
            Ap.x = A.xx*p.x + A.xy*p.y + A.xz*p.z;
            Ap.y = A.yx*p.x + A.yy*p.y + A.yz*p.z;
            Ap.z = A.zx*p.x + A.zy*p.y + A.zz*p.z;
            r.x = r0.x + Ap.x; r.y = r0.y + Ap.y; r.z = r0.z + Ap.z;
            double r2 = r.x*r.x + r.y*r.y + r.z*r.z;
            if (r2 < a_0*a_0){
              pcount++;
              double s = r2/(a_0*a_0) - 1.;
              // add V(|r|) = 1/|r| - tau(s)/a_0 to sum
              // F on particle i = V'(|r|)r/|r|
              // = (- 1/|r|^3 - 2 tau'(s)/a_0^3)r
              double taus = ff->tau[nu-1];
              double tau1s = 0.;
              for (int k = nu - 2; k >= 0; k--){
                tau1s = taus + s*tau1s;
                taus = ff->tau[k] + s*taus;}
              sum += 1./sqrt(r2) - taus/a_0;
              double coeff = -1./(r2*sqrt(r2)) - 2.*tau1s/pow(a_0, 3);
              F.x += coeff*r.x, F.y += coeff*r.y, F.z += coeff*r.z;}}
      
      energy += charge[i]*charge[j]*sum;
      F.x *= charge[i]*charge[j];
      F.y *= charge[i]*charge[j];
      F.z *= charge[i]*charge[j];
      force[i].x += F.x, force[i].y += F.y,force[i].z += F.z;
      force[j].x -= F.x, force[j].y -= F.y,force[j].z -= F.z;


      j = nlist[ptr++];
      if (pcount > 0) int_list++;
    }
    i = nlist[ptr++];
  }
  //printf("int_list : %d\n",int_list);
  return energy;}

static double Bspline(FF *ff, int i, double t);
static void anterpolate(FF *ff, Triple gd, double *q, int N, double *charge,
                        Vector *r){

  Matrix Ai = *(Matrix *)ff->Ai;
  Triple sd = {0,0,0};
  Triple of = {0,0,0};
  int nu = ff->orderAcc;
  if (ff->maxLevel > 1) { // q has padding
    int dmax = 2*ff->nLim + 1;
    sd.x = max(min(dmax,gd.x),nu-1);
    sd.y = max(min(dmax,gd.y),nu-1);
    sd.z = max(min(dmax,gd.z),nu-1);
    of.x = max(sd.x/2, nu/2 - 1);
    of.y = max(sd.y/2, nu/2 - 1);
    of.z = max(sd.z/2, nu/2 - 1);
  } else if (ff->maxLevel == 1) {
    sd.x = nu - 1;
    sd.y = nu - 1;
    sd.z = nu - 1;
    of.x = nu/2 - 1;
    of.y = nu/2 - 1;
    of.z = nu/2 - 1;
  }
  
  for (int i = 0; i < N; i++){
    Vector ri = r[i];
    Vector s = prod(Ai, ri);
    s.x = s.x - floor(s.x);
    s.y = s.y - floor(s.y);
    s.z = s.z - floor(s.z);
    Vector t = {(double)gd.x*s.x, (double)gd.y*s.y, (double)gd.z*s.z};
    Vector em = {floor(t.x), floor(t.y), floor(t.z)};
    t.x -= em.x, t.y -= em.y, t.z -= em.z;
    Triple m = {(int)em.x, (int)em.y, (int)em.z};
      for (int nx = - nu/2; nx < nu/2; nx++){
      double Qi = Bspline(ff, nx + nu/2, t.x);
        for (int ny = - nu/2; ny < nu/2; ny++){
        double Qj = Bspline(ff, ny + nu/2, t.y);
          for (int nz = - nu/2; nz < nu/2; nz++){
          double Qk = Bspline(ff, nz + nu/2, t.z);
          // add to q_{m+n}
          Triple m_n = {((m.x - nx) % gd.x + gd.x) % gd.x,
                        ((m.y - ny) % gd.y + gd.y) % gd.y,
                        ((m.z - nz) % gd.z + gd.z) % gd.z};
          q[((m_n.x + of.x)*(gd.y + sd.y) + (m_n.y + of.y))*(gd.z + sd.z) + m_n.z + of.z] += charge[i]*Qi*Qj*Qk;}}}}

  int gdznew = gd.z + sd.z;
  int gdynew = gd.y + sd.y;
  int gdxnew = gd.x + sd.x;
  for (int i = 0 ; i < gdxnew ;i++) {
    int qx = i;
    while (qx < of.x)            qx += gd.x;
    while (qx > gd.x - 1 + of.x) qx -= gd.x;
    for (int j = 0 ; j < gdynew ;j++) {
      int qy = j;
      while (qy < of.y)            qy += gd.y;
      while (qy > gd.y - 1 + of.y) qy -= gd.y;
      for (int k = 0 ; k < gdznew ;k++) {
        int qz = k;
        while (qz < of.z)            qz += gd.z;
        while (qz > gd.z - 1 + of.z) qz -= gd.z;
        
        // i,j,k -> outside
        if (i < of.x || i > of.x + gd.x - 1 ||
            j < of.y || j > of.y + gd.y - 1 ||
            k < of.z || k > of.z + gd.z - 1 ) {
          int inside  = qx * gdynew * gdznew + qy * gdznew + qz;
          int outside =  i * gdynew * gdznew +  j * gdznew + k;
          
          q[inside] += q[outside];
        }
      }
    }
  }
}


static void restrict_(FF *ff, Triple gd, double *ql, double *qlm1, int lm1){
  // gd are grid dimensions for ql, ql is initially zero
  // :::this can be made more efficient:::
  int dmax = 2*ff->nLim + 1;
  int nu = ff->orderAcc;
  Triple sd = {0,0,0};
  Triple of = {0,0,0};
  if (ff->maxLevel != lm1 + 1 ) { // ql has padding
    sd.x = max(min(dmax,gd.x),nu);
    sd.y = max(min(dmax,gd.y),nu);
    sd.z = max(min(dmax,gd.z),nu);
    of.x = max(sd.x/2, nu/2);
    of.y = max(sd.y/2, nu/2);
    of.z = max(sd.z/2, nu/2);
  }
  
  int tsx = max(min(dmax,2*gd.x),nu);
  int tsy = max(min(dmax,2*gd.y),nu);
  int tsz = max(min(dmax,2*gd.z),nu);
  int tdx = 2*gd.x + tsx;
  int tdy = 2*gd.y + tsy;
  int tdz = 2*gd.z + tsz;
  int tox = max(tsx/2, nu/2);
  int toy = max(tsy/2, nu/2);
  int toz = max(tsz/2, nu/2);
  double *J = ff->J + nu/2;

  /*
  for (int mx = 0; mx < gd.x; mx++)
    for (int my = 0; my < gd.y; my++)
      for (int mz = 0; mz < gd.z; mz++){
        int m = ((mx+of.x)*(gd.y + sd.y) + (my + of.y))*(gd.z + sd.z) + mz + of.z;
        for (int nx = - nu/2; nx <= nu/2; nx++){
          int tmpnx = 2*mx + nx + tox;
          for (int ny = - nu/2; ny <= nu/2; ny++){
            int tmpny = 2*my + ny + toy;
            for (int nz = - nu/2; nz <= nu/2; nz++){
              int tmpnz = 2*mz + nz + toz;
              int tmpn = (tmpnx*tdy + tmpny)*tdz + tmpnz;
              ql[m] += J[nx]*J[ny]*J[nz]*qlm1[tmpn];}}}} */

  double *F = (double *)calloc(tdx*tdy*(gd.z + sd.z),sizeof(double));
  double *G = (double *)calloc(tdx*(gd.y+sd.y)*(gd.z + sd.z),sizeof(double));
  // operator action on z-dimension
  for (int i = 0; i < tdx ; i++) {
    for (int j = 0 ; j < tdy ; j++) {
      for (int k = 0 ; k < gd.z ; k++) {
        double sum = 0.0;
        for (int nz = - nu/2; nz <= nu/2; nz++) {
          double value = qlm1[(i*tdy + j)*tdz + 2*k+nz + toz];
          sum += J[nz] * value; }
        F[(i*tdy + j)*(gd.z + sd.z)  + k + of.z] = sum; }}}

  // operator action on y-dimension
  for (int i = 0; i < tdx ; i++) {
    for (int j = 0 ; j < gd.y ; j++) {
      for (int k = 0 ; k < gd.z ; k++) {
        double sum = 0.0;
        for (int ny = - nu/2; ny <= nu/2; ny++) {
          double value = F[(i*tdy + 2*j+ny+toy)*(gd.z + sd.z)  + k + of.z] ;
          sum += J[ny] * value; }
        G[(i*(gd.y + sd.y) + j+of.y)*(gd.z + sd.z)  + k + of.z] = sum; }}}

  // operator action on x-dimension
   for (int i = 0; i < gd.x ; i++) {
    for (int j = 0 ; j < gd.y ; j++) {
      for (int k = 0 ; k < gd.z ; k++) {
        double sum = 0.0;
        for (int nx = - nu/2; nx <= nu/2; nx++) {
          double value = G[((2*i+nx+tox)*(gd.y + sd.y) + j+of.y)*(gd.z + sd.z)  + k + of.z]; 
          sum += J[nx] * value; }
        ql[((i+of.x)*(gd.y + sd.y) + j+of.y)*(gd.z + sd.z)  + k + of.z] = sum; }}}

  free(F);
  free(G);

}



static void DFT(Triple gd, double complex dL[], double fL[]);
static void invDFT(Triple gd, double fL[], double complex dL[]);
static void DFTg2g(Triple gd, double *el, double *ql, double *kh){
  // el = DFT of ql
  double complex *dl
    = (double complex *)malloc(gd.x*gd.y*gd.z*sizeof(double complex));
  DFT(gd, dl, ql);
  // ql = kh . el pointwise
  for (int i = 0; i < gd.x; i++)
    for (int j = 0; j < gd.y; j++)
      for (int k = 0; k < gd.z; k++){
        dl[(i*gd.y + j)*gd.z + k] *= kh[(i*gd.y + j)*gd.z + k];
      }
  // el = invDFT of ql
  invDFT(gd, el, dl);
  free(dl);
}

static void FFTg2g(FF *ff, Triple gd, double *el, double *ql, double *kh){
  int nu = ff->orderAcc;
#ifdef NO_FFT
  ;
#else
  int ofx = nu/2 - 1;
  int ofy = nu/2 - 1;
  int ofz = nu/2 - 1;
  
  // el = DFT of ql
  int c = 0;
  for (int i = 0; i < gd.x; i++) {
    for (int j = 0 ; j < gd.y; j++) {
      for (int k = 0 ; k < gd.z ; k++) {
        ff->fftw_in[c++] = (fftw_complex)
           ql[((i+ofx)*(gd.y + nu-1) + (j+ofy))*(gd.z + nu - 1) + k + ofz];
      }
    }
  }
  
  fftw_execute(ff->forward);
  // ql = kh . el pointwise
  for (int i = 0; i < gd.x*gd.y*gd.z; i++)
        ff->fftw_in[i] *= (fftw_complex)kh[i];
  // el = invDFT of ql
  fftw_execute(ff->backward);
  for (int i = 0; i < gd.x*gd.y*gd.z; i++)
    el[i] = creal(ff->fftw_in[i])/(double)(gd.x*gd.y*gd.z);
#endif
}

static void prolongate(FF *ff, Triple gd, double *el, double *elp1){
  // add J_n elp1_m to el_{2m+n}
  int nu = ff->orderAcc;
  int hdx = gd.x/2, hdy = gd.y/2, hdz = gd.z/2;
  double *J = ff->J + nu/2;

  double *F = (double *)calloc(hdx*hdy*gd.z, sizeof(double));
  double *G = (double *)calloc(hdx*gd.y*gd.z, sizeof(double));
  /* Along z-dimension */
  for (int mx = 0; mx < hdx; mx++){
    for (int my = 0; my < hdy; my++){
      for (int mz = 0; mz < hdz; mz++){
        double value = elp1[(mx*hdy + my)*hdz + mz];
        for (int nz = - nu/2; nz <= nu/2; nz++){
          int kz = ((2*mz + nz) % gd.z + gd.z) % gd.z;
          F[(mx*hdy + my)*gd.z + kz] += J[nz]*value;}}}}

  /* Along y-dimension */
   for (int mx = 0; mx < hdx; mx++){
    for (int my = 0; my < hdy; my++){
      for (int mz = 0; mz < gd.z; mz++){
        double value = F[(mx*hdy + my)*gd.z + mz];
        for (int ny = - nu/2; ny <= nu/2; ny++){
          int ky = ((2*my + ny) % gd.y + gd.y) % gd.y;
          G[(mx*gd.y + ky)*gd.z + mz] += J[ny]*value;}}}}

   /* Along z-dimension */
   for (int mx = 0; mx < hdx; mx++){
    for (int my = 0; my < gd.y; my++){
      for (int mz = 0; mz < gd.z; mz++){
        double value = G[(mx*gd.y+ my)*gd.z + mz];
        for (int nx = - nu/2; nx <= nu/2; nx++){
          int kx = ((2*mx + nx) % gd.x + gd.x) % gd.x;
          el[(kx*gd.y + my)*gd.z + mz] += J[nx]*value;}}}}
  
  free(F);
  free(G);
}

static double Bspline1(double *Qi1t, FF *ff, int i, double t);
static void interpolate(FF *ff, int N, Vector *E, Vector *r, Triple gd,
 double *el){
  // return grid-level electric field
  Matrix Ai = *(Matrix *)ff->Ai;
  int nu = ff->orderAcc;
  for (int i = 0; i < N; i++){
    double Ex = 0., Ey = 0., Ez = 0.;
    Vector ri = r[i];
    Vector s = prod(Ai, ri);
    Vector t = {(double)gd.x*s.x, (double)gd.y*s.y, (double)gd.z*s.z};
    Vector em = {floor(t.x), floor(t.y), floor(t.z)};
    t.x -= em.x, t.y -= em.y, t.z -= em.z;
    Triple m = {(int)em.x, (int)em.y, (int)em.z};
    //**for (int ni = 0; ni < nu; ni++){
    for (int ni = - nu/2; ni < nu/2; ni++){
      double Qi1;
      //**double Qi = Bspline1(&Qi1, ff, ni, t.x);
      double Qi = Bspline1(&Qi1, ff, ni + nu/2, t.x);
      //**for (int nj = 0; nj < nu; nj++){
      for (int nj = - nu/2; nj < nu/2; nj++){
        double Qj1;
        //**double Qj = Bspline1(&Qj1, ff, nj, t.y);
        double Qj = Bspline1(&Qj1, ff, nj + nu/2, t.y);
        //**for (int nk = 0; nk < nu; nk++){
        for (int nk = - nu/2; nk < nu/2; nk++){
          double Qk1;
          //**double Qk = Bspline1(&Qk1, ff, nk, t.z);
          double Qk = Bspline1(&Qk1, ff, nk + nu/2, t.z);
          Triple m_n = {((m.x - ni) % gd.x + gd.x) % gd.x,
                        ((m.y - nj) % gd.y + gd.y) % gd.y,
                        ((m.z - nk) % gd.z + gd.z) % gd.z};
          double elm_n = el[(m_n.x*gd.y + m_n.y)*gd.z + m_n.z];
          Ex -= Qi1*Qj*Qk*elm_n;
          Ey -= Qi*Qj1*Qk*elm_n;
          Ez -= Qi*Qj*Qk1*elm_n;}}}
    Ex *= (double)gd.x, Ey *= (double)gd.y, Ez *= (double)gd.z;
    Vector E_ = {Ex, Ey, Ez};
    E[i] = prodT(Ai, E_);}}

static double Bspline(FF *ff, int i, double t){
  // assumes 0 <= i < nu and 0 <= t < 1
  int nu = ff->orderAcc;
  if (i >= nu/2){i = nu - 1 - i; t = 1. - t;}
  // piece_i(i + t) = q_i0 + q_i1*t + ... + q_{i,nu-1}*t^{nu-1}
  double *Qi = ff->Q + i*nu;
  double Qit = Qi[nu-1];
  for (int k = nu - 2; k >= 0; k--)
    Qit = Qi[k] + t*Qit;
  //static int counter = 0;
  //if (i >1 ){
  //  printf("nu:%d i:%d t:%f Qit:%f\n",nu,i,t,Qit);
  //  counter++;}
  return Qit;}

static double Bspline1(double *Qi1t, FF *ff, int i, double t){
  // assumes 0 <= i < nu and 0 <= t < 1
  // piece_i(i + t) = q_i0 + q_i1*t + ... + q_{i,nu-1}*t^{nu-1}
  int nu = ff->orderAcc;
  double sign = 1.;
  if (i >= nu/2){i = nu - 1 - i; t = 1. - t; sign = -1.;}
  double *Qi = ff->Q + i*nu;
  double Qit = Qi[nu-1];
  *Qi1t = 0.;
  for (int k = nu - 2; k >= 0; k--){
    *Qi1t = Qit + t*(*Qi1t);
    Qit = Qi[k] + t*Qit;}
  *Qi1t *= sign;
  return Qit;}

Vector prod(Matrix M, Vector v){
  Vector Mx = {M.xx*v.x + M.xy*v.y + M.xz*v.z,
               M.yx*v.x + M.yy*v.y + M.yz*v.z,
               M.zx*v.x + M.zy*v.y + M.zz*v.z};
  return Mx;}

Vector prodT(Matrix M, Vector v){
  Vector MTx = {M.xx*v.x + M.yx*v.y + M.zx*v.z,
                M.xy*v.x + M.yy*v.y + M.zy*v.z,
                M.xz*v.x + M.yz*v.y + M.zz*v.z};
  return MTx;}

static void DFT(Triple gd, double complex dL[], double fL[]){
  double twopi = 8.*atan(1.);
  for (int kx = 0; kx < gd.x; kx++)
    for (int ky = 0; ky < gd.y; ky++)
      for (int kz = 0; kz < gd.z; kz++){
        double complex dLk = 0.;
        double cx = twopi*(double)kx/gd.x;
        double cy = twopi*(double)ky/gd.y;
        double cz = twopi*(double)kz/gd.z;
        for (int nx = 0; nx < gd.x; nx++)
          for (int ny = 0.; ny < gd.y; ny++)
            for (int nz = 0.; nz < gd.z; nz++){
              double cn = cx*(double)nx + cy*(double)ny + cz*(double)nz;
              dLk += cexp(-cn*I)*fL[(nx*gd.y + ny)*gd.z + nz];}
        dL[(kx*gd.y + ky)*gd.z + kz] = dLk;
      }
}

static void invDFT(Triple gd, double fL[], double complex dL[]){
  double twopi = 8.*atan(1.);
  for (int nx = 0; nx < gd.x; nx++)
    for (int ny = 0; ny < gd.y; ny++)
      for (int nz = 0; nz < gd.z; nz++){
        double fLn = 0.;
        double cx = twopi*(double)nx/gd.x;
        double cy = twopi*(double)ny/gd.y;
        double cz = twopi*(double)nz/gd.z;
        for (int kx = 0; kx < gd.x; kx++)
          for (int ky = 0.; ky < gd.y; ky++)
            for (int kz = 0.; kz < gd.z; kz++){
              double ck = cx*(double)kx + cy*(double)ky + cz*(double)kz;
              fLn += creal(cexp(ck*I)*dL[(kx*gd.y + ky)*gd.z + kz]);}
        fL[(nx*gd.y + ny)*gd.z + nz] = fLn/(gd.x*gd.y*gd.z);
      }
}

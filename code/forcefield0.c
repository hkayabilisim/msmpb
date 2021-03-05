// file forcefield.c
// compile with -lfftw3 or -DNO_FFT

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
#include <time.h>
#include "fftw3.h"
#include "forcefield.h"


void neighborlist(FF *ff, int N, Vector *position){
  int interaction_count = 0;
  Matrix A = *(Matrix *)ff->A;
  Matrix Ai = *(Matrix *)ff->Ai;
  double a_0 = ff->aCut[0];
  double a_02 = a_0 * a_0;
  Vector as = {sqrt(Ai.xx*Ai.xx + Ai.yx*Ai.yx + Ai.zx*Ai.zx),
    sqrt(Ai.xy*Ai.xy + Ai.yy*Ai.yy + Ai.zy*Ai.zy),
    sqrt(Ai.xz*Ai.xz + Ai.yz*Ai.yz + Ai.zz*Ai.zz)};
  double pi = 4.*atan(1.);
  int nlist_idx = 0;
  int nlist_len = 3*N + 2.1 * ((4.0/3.0)*pi*a_0*a_0*a_0*0.5*N*(N-1.0)/ff->detA);
  ff->nlist_len = nlist_len;
  int *nlist = (int *) calloc(nlist_len, sizeof(int));
  ff->nlist = nlist;
  ff->nlist_len = nlist_len;
  
  // create hash table
  Triple gd = {ceil(pow(N,1./3.)),ceil(pow(N,1./3.)),ceil(pow(N,1./3.))};
  int gd_prod = gd.x*gd.y*gd.z;
  int *first = (int *)malloc(gd_prod*sizeof(int));
  for (int m = 0; m < gd_prod; m++) first[m] = -1;
  int *next = (int *)malloc(N*sizeof(int));
  for (int i = 0; i < N; i++){
    Vector ri = position[i];
    Vector s = prod(Ai, ri);
    s.x = s.x - floor(s.x);
    s.y = s.y - floor(s.y);
    s.z = s.z - floor(s.z);
    Vector t = {(double)gd.x*s.x, (double)gd.y*s.y, (double)gd.z*s.z};
    // Triple m = {(int)floor(t.x), (int)floor(t.y), (int)floor(t.z)};
    Triple m = {(int)t.x, (int)t.y, (int)t.z};
    int m_ = (m.x*gd.y + m.y)*gd.z + m.z;
    next[i] = first[m_];
    first[m_] = i;
  }
  // compute ranges - can be part of preprocessing
  // compute diameter of a grid cell
  Vector HAx = {A.xx/(double)gd.x, A.yx/(double)gd.y, A.zx/(double)gd.z};
  Vector HAy = {A.xy/(double)gd.x, A.yy/(double)gd.y, A.zy/(double)gd.z};
  Vector HAz = {A.xz/(double)gd.x, A.yz/(double)gd.y, A.zz/(double)gd.z};
  double off_yz = HAy.x*HAz.x + HAy.y*HAz.y + HAy.z*HAz.z;
  double off_zx = HAz.x*HAx.x + HAz.y*HAx.y + HAz.z*HAx.z;
  double off_xy = HAx.x*HAy.x + HAx.y*HAy.y + HAx.z*HAy.z;
  if (off_yz * off_zx * off_xy < 0){
    if (fabs(off_yz) < fabs(off_zx) && fabs(off_yz) < fabs(off_xy))
      off_yz *= -1.;
    else if (fabs(off_zx) < fabs(off_xy))
      off_zx *= -1.;
    else
      off_xy *= -1.;}
  Vector t = {copysign(1., off_yz), copysign(1., off_zx),copysign(1., off_xy)};
  Vector At = prod(A, t);
  Vector HAt = {At.x/(double)gd.x, At.y/(double)gd.y, At.z/(double)gd.z};
  double diam = sqrt(HAt.x*HAt.x + HAt.y*HAt.y + HAt.z*HAt.z);
  int nxlim = (int)ceil(as.x*(double)gd.x*a_0);
  int nylim = (int)ceil(as.y*(double)gd.y*a_0);
  int nzlim = (int)ceil(as.z*(double)gd.z*a_0);
  
  int nxd = 2*nxlim + 1 < gd.x ? 2*nxlim + 1 : gd.x;
  int nyd = 2*nylim + 1 < gd.y ? 2*nylim + 1 : gd.y;
  int nzd = 2*nzlim + 1 < gd.z ? 2*nzlim + 1 : gd.z;
  // compute cell neighbor lists - can be part of preprocessing
  int ndim = nxd*nyd*nzd;
  Triple *n = (Triple *)calloc(ndim, sizeof(Triple));
  int k = 0; // number of occupied elements of n;
  for (int nx = -nxd/2; nx < (nxd + 1)/2; nx++)
    for (int ny = -nyd/2; ny < (nyd + 1)/2; ny++)
      for (int nz = -nzd/2; nz < (nzd + 1)/2; nz++){
        Vector s = {(double)nx/(double)gd.x,
          (double)ny/(double)gd.y, (double)nz/(double)gd.z};
        Vector r = prod(A, s);
        double normr = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
        if (normr - diam < a_0){
          Triple nk = {nx, ny, nz};
          n[k] = nk;
          k++;}
      }
  int kcnt = k;
  // compute interactions
  int m = 0;  // index of grid cell 1
  int i = -1;  // particle number
  while (true){  // loop over i
    if (i >= 0) i = next[i];
    else i = first[m];
    while (i < 0 && m + 1 < gd_prod){m++; i = first[m];}
    //### loop exit ###
    if (i < 0) break;
    int mx = m/(gd.y*gd.z);
    int my = (m - mx*gd.y*gd.z)/gd.z;
    int mz = m - (mx*gd.y + my)*gd.z;
    Vector ri = position[i];
    nlist[nlist_idx++] = i;
    int k = 0; // index of next neighbor triple for j
    int j = i;
    while (true){  // loop over j
      j = next[j];
      while (j < 0 && k < kcnt){
        int mpnx = (mx + n[k].x + gd.x)%gd.x;
        int mpny = (my + n[k].y + gd.y)%gd.y;
        int mpnz = (mz + n[k].z + gd.z)%gd.z;
        int mpn = (mpnx*gd.y + mpny)*gd.z + mpnz;
        if (mpn > m) j = first[mpn];
        k++;}
      //### loop exit ###
      if (j < 0) break;
      // compute interaction
      Vector rj = position[j];
      Vector r = {rj.x - ri.x, rj.y - ri.y, rj.z - ri.z};
      // convert to nearest image
      Vector s = prod(Ai, r);
      Vector p = {floor(s.x + 0.5), floor(s.y + 0.5), floor(s.z + 0.5)};
      Vector Ap = prod(A, p);
      r.x -= Ap.x; r.y -= Ap.y; r.z -= Ap.z;
      double distance2 = r.x*r.x + r.y*r.y + r.z*r.z;
      if (distance2 < a_02) {
        interaction_count++;
        nlist[nlist_idx++] = j;
      }
    }
    nlist[nlist_idx++] = -1;
  }
  nlist[nlist_idx++] = -1;
  free(n);
  free(first);
  free(next);
  ff->nlist_interaction_count = interaction_count;
}

FF *FF_new(int N, double *q, double edges[3][3]){
  FF *ff = (FF *)calloc(1, sizeof(FF));
#ifndef NO_FFT
  ff->FFT = true;
#endif
  ff->N = N;
  ff->q = q;
  for (int i = 0 ; i < 3; i++)
    for (int j = 0 ; j < 3; j++)
      ff->A[i][j] = edges[i][j];
  ff->maxLevel = -1;
  ff->orderAcc = -1;
  ff->topGridDim[0] = -1;
  ff->topGridDim[1] = -1;
  ff->topGridDim[2] = -1;
  ff->cutoff = -1;
  ff->estErr = -1;
  return ff;}

// for each computational parameter
// have a get method FF_set_... that returns the parameter
// Each set method is optional
// No method parameter should be changed after a call to build
// ::: YET TO DO: PARAMETER CHECKING :::
void FF_set_orderAcc(FF *ff, int orderAcc ){
  // orderAcc: even integer >= 4
  ff->orderAcc = orderAcc;}
void FF_set_maxLevel(FF *ff, int maxLevel) {
  // level-l cutoff relative to level-l grid size
  //: currently relative to minimum box dimension
  ff->maxLevel = maxLevel;}
void FF_set_topGridDim(FF *ff, int topGridDim[3]){
  *(Triple *)ff->topGridDim = *(Triple *)topGridDim;}

void FF_set_estErr(FF *ff,double estErr){
  ff->estErr = estErr;
}
double FF_get_estErr(FF *ff){
  return ff->estErr;
}
void FF_set_cutoff(FF *ff,double cutoff) {
  ff->cutoff = cutoff;
}
double FF_get_cutoff(FF *ff) {
  return ff->cutoff;
}
void FF_set_FFT(FF *ff, bool FFT){
#ifdef NO_FFT
  if (FFT){
    printf("FFT unavailable; execution aborted\n");
    exit(1);}
#endif
  ff->FFT = FFT;}

// helper functions:

static int increment(int Minput){
  int M = Minput;
  do {
    int factor = 2;
    int remainder = M;
    bool hasnonprimefactor = false;
    do {
      if (remainder % factor == 0) {
        remainder = remainder / factor;
        if (factor != 2 && factor != 3 && factor != 5 && factor != 7) {
          hasnonprimefactor = true;
          break;
        }
      } else
        factor++;
    } while (remainder > 1);
    
    if (hasnonprimefactor)
      M++;
    else {
      //if (Minput != M) {
        //fprintf(stderr,"warning M changed from %d to %d\n",Minput,M);
      //}
      return M;
    }
  } while (1);
}

void determineM(FF *ff,int L,int *Mout) {
  Matrix A = *(Matrix *)ff->A;
  int N = ff->N;
  double eta = 1.0;
  if (ff->FFT) {
    if (L >=2)
      eta = 0.83;
    else
      eta = 5.7 ;
  }

  double ax = sqrt(A.xx*A.xx + A.xy*A.xy + A.xz*A.xz);
  double ay = sqrt(A.yx*A.yx + A.yy*A.yy + A.yz*A.yz);
  double az = sqrt(A.zx*A.zx + A.zy*A.zy + A.zz*A.zz);
  
  double Msx = pow(2.,1.-L)*pow(eta*N,1./3.)*ax*pow(ax*ay*az,-1./3.);
  double Msy = pow(2.,1.-L)*pow(eta*N,1./3.)*ay*pow(ax*ay*az,-1./3.);
  double Msz = pow(2.,1.-L)*pow(eta*N,1./3.)*az*pow(ax*ay*az,-1./3.);
  
  Mout[0] = ceil(Msx);
  Mout[1] = ceil(Msy);
  Mout[2] = ceil(Msz);
}

void FF_build(FF *ff, double (*position)[3]){
  // N: number of particles
  // edges: columns are vectors defining parallelpiped
  int N = ff->N;
  for (int i = 0 ; i < 3; i++)
    for (int j = 0 ; j < 3; j++)
      ff->Ai[i][j] = ff->A[i][j];
  Matrix A = *(Matrix *)ff->A;
  Matrix Ai = *(Matrix *)ff->Ai;

  ff->detA = invert(&Ai);
  *(Matrix *)ff->Ai = Ai;
  
  // Default error tolerance
  if (ff->estErr == -1) ff->estErr = 0.001;
  
  // Set L=1 if not specified
  if (ff->maxLevel == -1) ff->maxLevel = 1;
  
  if (! ff->FFT &&  ff->topGridDim[0] == -1) {
    int l = 1 ;
    int Mest[3];
    do {
      determineM(ff,l,&Mest[0]);
      int Mprod = Mest[0]*Mest[1]*Mest[2];
      if (Mprod <= sqrt(N))  break;
      l++;
    } while (true);
    if (ff->maxLevel != l)
      fprintf(stderr,"L=%d is replaced by L=%d\n",ff->maxLevel,l);
    ff->maxLevel = l;
    ff->topGridDim[0] = Mest[0];
    ff->topGridDim[1] = Mest[1];
    ff->topGridDim[2] = Mest[2];
  }
  
  if (ff->FFT && ff->topGridDim[0] == -1) {
    determineM(ff,ff->maxLevel,ff->topGridDim);
  }
  
  if (ff->FFT) {
    ff->topGridDim[0] = increment(ff->topGridDim[0]);
    ff->topGridDim[1] = increment(ff->topGridDim[1]);
    ff->topGridDim[2] = increment(ff->topGridDim[2]);
  }
  
  // determine a0 if not specified
  if (ff->cutoff == -1) {
    if (ff->orderAcc == -1) { // order not selected
      double min_cutoff = INFINITY;
      for (int nu_test = 4 ; nu_test <= 10; nu_test += 2) {
        double top = _FF_get_errEst(ff,nu_test);
        double cutoff_test = pow(top/ff->estErr,1./(nu_test));
        if (cutoff_test < min_cutoff)
          min_cutoff = cutoff_test;
      }
      ff->cutoff = min_cutoff;
    } else {
      double top = _FF_get_errEst(ff,ff->orderAcc);
      ff->cutoff = pow(top/ff->estErr,1./(ff->orderAcc));
    }
  }
  
  // determine nu
  if (ff->orderAcc == -1) {
    double min_err = INFINITY;
    int best_nu = 4;
    for (int nu_test = 4 ; nu_test <= 10; nu_test += 2) {
      double err = _FF_get_errEst(ff,nu_test)/pow(ff->cutoff,nu_test);
      if (err < min_err) {
        min_err = err;
        best_nu = nu_test;
      }
    }
    ff->orderAcc = best_nu;
  }
  




  ff->aCut = (double *)malloc((ff->maxLevel + 1)*sizeof(double));
  ff->aCut[0] = ff->cutoff;
  for (int l = 1; l <= ff->maxLevel; l++)
    ff->aCut[l] = 2.0*ff->aCut[l-1];
  
  // calculate relative cutoff
  double asx = sqrt(Ai.xx*Ai.xx + Ai.xy*Ai.xy + Ai.xz*Ai.xz);
  double asy = sqrt(Ai.yx*Ai.yx + Ai.yy*Ai.yy + Ai.yz*Ai.yz);
  double asz = sqrt(Ai.zx*Ai.zx + Ai.zy*Ai.zy + Ai.zz*Ai.zz);
  int nbarx = ceil(ff->aCut[ff->maxLevel] * ff->topGridDim[0] * asx);
  int nbary = ceil(ff->aCut[ff->maxLevel] * ff->topGridDim[1] * asy);
  int nbarz = ceil(ff->aCut[ff->maxLevel] * ff->topGridDim[2] * asz);
  ff->relCutoff = fmax(fmax(nbarx,nbary),nbarz);
  ff->nLim = ceil(ff->relCutoff - 1.);

  // build tau(s) = tau_0 + tau_1*s + ... + tau_{ord-1} *s^{ord-1}
  ff->tau = (double *)malloc(ff->orderAcc*sizeof(double));
  ff->tau[0] = 1.;
  for (int i = 1; i < ff->orderAcc; i++)
    ff->tau[i] = - (1. - 0.5/(double)i)*ff->tau[i-1];
  
  // build first ord/2 pieces of B-splines Q_ord(t)
  // piece_i(t) = q_i0 + q_i1*(t - i) + ... + q_{i,ord-1}*(t - i)^{ord-1}
  int nu = ff->orderAcc, nknots = nu/2;
  ff->Q = (double *)calloc(nknots*nu, sizeof(double));
  // for degrees i = 0, 1, ..., nu-1
  // compute values of Q_{i+1} at knots j = 0, 1, ..., nknots-1
  double *Q = (double *)calloc(nu*nknots, sizeof(double));
  // Q_1(0) = 1, Q_1(j) = 0 for j > 0
  Q[0] = 1.;
  double *Qip1 = Q;
  // Q_{i+1}(j) = (j*Q_i(j) + (i+1-j)*Q_i(j-1))/i
  for (int i = 1; i < nu; i++){
    double *Qi = Qip1;
    Qip1 = Qi + nknots;
    Qip1[0] = 0.;
    for (int j = 1; j < nknots; j++)
      Qip1[j] =
        ((double)j*Qi[j] + (double)(i + 1 - j)*Qi[j-1])/(double)i;}
  // compute k-th derivative of Q_{k+1}, ..., Q_nu divided by k!
  // Q_{i+1}^(k)(j+)/k! = (Q_i^(k-1)(j+)/(k-1)! - Q_i^(k-1)(j-1+)/(k-1)!)/k
  for (int k = 1; k < nu; k++){
    // (k-1)st derivative of Q_{i+k} is in row i
    // replace with k-th derivative of Q{i+k+1}
    for (int i = 0; i < nu - k; i++){
      double *Qi = Q + i*nknots;
      for (int j = nknots - 1; j > 0; j--)
        Qi[j] = (Qi[j] - Qi[j-1])/(double)k;
      Qi[0] = Qi[0]/(double)k;}}
  for (int j = 0; j < nknots; j++)
    for (int k = 0; k < nu; k++)
      ff->Q[j*nu + k] = Q[(nu - 1 - k)*nknots + j];
  free(Q);

  // build two-scale stencil J[n], n = -nu/2, ..., nu/2
  ff->J = (double *)calloc(nu + 1, sizeof(double));
  // calculate Pascal's triangle
  double *J0 = ff->J + nu/2;
  J0[0] = 1.;
  for (int i = 1; i <= nu/2; i++){
    for (int j = -i; j < i; j++) J0[j] += J0[j+1];
    for (int j = i; j > -i; j--) J0[j] += J0[j-1];}
  for (int i = -nu/2; i <= nu/2; i++) J0[i] *= pow(2., 1 - nu);
  
  // d^i sigma (1-) = d^i gama(1-) + i * d^(i-1) gama(1-)
  // Implementation of A.1 from periodic.pdf (version 20170811)
  double *dsigma = malloc(sizeof(double) * (2*nu + 2));
  double **dgama = (double **)malloc(sizeof(double *)*(2*nu+2));
  for (int i=0 ; i<2*nu+2; i++)
    dgama[i] = (double *)malloc(sizeof(double)*nu);
  
  for (int k=nu-1;k>=0 ;k--) dgama[0][k] = 1;
  for (int i = 1; i < 2*nu ; i++) {
    dgama[i][nu-1] = 0;
    for (int k = nu - 1 ; k >= 1; k--) {
      dgama[i][k-1] = ((.5 - k)/k) * (2.*i*dgama[i-1][k] +
                                      (i > 1 ? (i*(i-1.0)*dgama[i-2][k]) : 0.0));
    }
  }
  dsigma[0] = dgama[0][0];
  for (int i = 1 ; i <= 2*nu + 1 ; i++)
    dsigma[i] = dgama[i][0] + i * dgama[i-1][0];
  ff->dsigma = dsigma;
  for (int i=0 ; i<2*nu+2; i++) free(dgama[i]);
  free(dgama);

  // build anti-blurring operator
  omegap(ff);
  
  // build grid-to-grid stencil
  // assume DFT used
  // allocate memory for khat
  // special treatment for level l = maxLevel
  int l = ff->maxLevel;
  ff->khat = (double **)malloc((l+1)*sizeof(double *));
  int di = ff->topGridDim[0], dj = ff->topGridDim[1], dk = ff->topGridDim[2];
  ff->khat[l] = (double *)calloc(di*dj*dk, sizeof(double));
#ifdef NO_FFT
  ff->fftw_in = (double complex *)malloc(di*dj*dk*sizeof(double complex));
#else
  ff->fftw_in = (fftw_complex *)fftw_malloc(di*dj*dk*sizeof(fftw_complex));
  ff->forward = fftw_plan_dft_3d(di, dj, dk, ff->fftw_in, ff->fftw_in,
                                FFTW_FORWARD, FFTW_ESTIMATE);
  ff->backward = fftw_plan_dft_3d(di, dj, dk, ff->fftw_in, ff->fftw_in,
                                FFTW_BACKWARD, FFTW_ESTIMATE);
#endif

  int dmax = 2*ff->nLim + 1;
  for (l = ff->maxLevel-1; l > 0; l--){
    di *= 2; dj *= 2; dk *= 2;
    int sdi = dmax < di ? dmax : di;
    int sdj = dmax < dj ? dmax : dj;
    int sdk = dmax < dk ? dmax : dk;
    ff->khat[l] = (double *)calloc(sdi*sdj*sdk, sizeof(double));}
  // following is for purpose of debugging
  if (strcmp(o.test, "nobuild") == 0) return;
  // compute stencils ff->khat[l]
  
  
  Vector *r = (Vector *)position;
  msm4g_tic();
  neighborlist(ff, N, r);
  ff->time_nlist = msm4g_toc();
  FF_rebuild(ff,ff->A,position);
}

// for each computational parameter
// also have a get method FF_set_... that returns the parameter
int FF_get_orderAcc(FF*ff) {
  return ff->orderAcc;}
int FF_get_maxLevel(FF *ff) {
  return ff->maxLevel;}
void FF_get_topGridDim(FF*ff, int topGridDim[3]) {
  *(Triple *)topGridDim = *(Triple *)ff->topGridDim;}
bool FF_get_FFT(FF *ff){
        return ff->FFT;}

double FF_get_deltaF(FF *ff) {
  // calculate C_{nu-1}
  int N = ff->N;
  double *charge = ff->q;
  int L = ff->maxLevel;
  int nu = ff->orderAcc;
  double a0 = ff->aCut[0];

  int index = (nu - 4)/2;
  //double C[4] = {0.21,1.24,17.59,452.94};
  double C[4] = {0.825169,5.16701,73.6508,1977.23};
  
  Matrix A = *(Matrix *)ff->A;
  double ax = sqrt(A.xx*A.xx + A.xy*A.xy + A.xz*A.xz);
  double ay = sqrt(A.yx*A.yx + A.yy*A.yy + A.yz*A.yz);
  double az = sqrt(A.zx*A.zx + A.zy*A.zy + A.zz*A.zz);
  Triple M1 = *(Triple *)ff->topGridDim;
  for (int i = 1; i < L; i++) {
    M1.x *= 2;
    M1.y *= 2;
    M1.z *= 2;
  }
  // complete the calculation
  double Q2 = 0;
  for (int i = 0; i < N; i++)
    Q2 += charge[i]*charge[i];
  
  double h1 = fmax(fmax(ax/M1.x,ay/M1.y),az/M1.z);
  double deltaF = C[index]*Q2*(1-pow(4.0,-L))*pow(h1,nu-1)*pow(ff->detA,-1./3.)*pow(N,-2.0/3.0)/pow(a0,nu);
  return deltaF;
}

double _FF_get_errEst(FF *ff,int nu) {
  // calculate C_{nu-1}
  int N = ff->N;
 double *charge = ff->q;
  int L = ff->maxLevel;

  int index = (nu - 4)/2;
  //double C[4] = {0.21,1.24,17.59,452.94};
  double C[4] = {1,1,1,1};
  
  Matrix A = *(Matrix *)ff->A;
  double ax = sqrt(A.xx*A.xx + A.xy*A.xy + A.xz*A.xz);
  double ay = sqrt(A.yx*A.yx + A.yy*A.yy + A.yz*A.yz);
  double az = sqrt(A.zx*A.zx + A.zy*A.zy + A.zz*A.zz);
  Triple M1 = *(Triple *)ff->topGridDim;
  for (int i = 1; i < L; i++) {
    M1.x *= 2;
    M1.y *= 2;
    M1.z *= 2;
  }
  // complete the calculation
  double Q2 = 0;
  for (int i = 0; i < N; i++)
    Q2 += charge[i]*charge[i];
  double Fref = (Q2/(double)N)/pow(ff->detA/(double)N, 5./6.);

  double h1 = fmax(fmax(ax/M1.x,ay/M1.y),az/M1.z);
  double deltaF = C[index]*Q2*(1-pow(4.0,-L))*pow(h1,nu-1)*pow(ff->detA,-1./3.)*pow(N,-2.0/3.0);
  return deltaF/Fref;
}
double FF_get_errEst(FF *ff){
  int nu = ff->orderAcc;
  double a0 = ff->aCut[0];
  return _FF_get_errEst(ff,nu)/pow(a0,nu);
}


// The rebuild method is called every time the periodic cell changes.
// It initializes edges and calculate the grid2grid stencils.
// helper functions:
static void dAL(FF *ff, Triple gd, Triple sd, double kh[], double detA);
static void kaphatA(FF *ff, int l, Triple gd, Triple sd, double kh[],
                    Vector as);
static void DFT(Triple gd, double dL[], double khatL[]);
static void FFT(FF *ff, Triple gd, double dL[], double khatL[]);
void FF_rebuild(FF *ff, double edges[3][3], double (*r)[3]) {
  *(Matrix *)ff->A = *(Matrix *)edges;
  Matrix A = *(Matrix *)ff->A;
  // calculate coeffs for const part of energy
  double *tau = ff->tau;
  int nu = ff->orderAcc;
  double a_0 = ff->aCut[0];
  int L = ff->maxLevel;
  double aL=ff->aCut[L];
  Matrix Ai = *(Matrix *)edges;
  double detA = invert(&Ai);
  *(Matrix *)ff->Ai = Ai;
  ff->detA = detA ;
  Vector as = {sqrt(Ai.xx*Ai.xx + Ai.yx*Ai.yx + Ai.zx*Ai.zx),
               sqrt(Ai.xy*Ai.xy + Ai.yy*Ai.yy + Ai.zy*Ai.zy),
               sqrt(Ai.xz*Ai.xz + Ai.yz*Ai.yz + Ai.zz*Ai.zz)};
  double xlim = ceil(a_0*as.x - 1.);
  double ylim = ceil(a_0*as.y - 1.);
  double zlim = ceil(a_0*as.z - 1.);
  ff->coeff1 = 0.;  // coeff of 1/2(sum_i q_i^2)
  for (double x = -xlim; x <= xlim; x++)
    for (double y = -ylim; y <= ylim; y++)
      for (double z = -zlim; z <= zlim; z++){ // could do just z >= 0
        if (x*x + y*y + z*z == 0) continue;
        double rx = A.xx*x + A.xy*y + A.xz*z;
        double ry = A.yx*x + A.yy*y + A.yz*z;
        double rz = A.zx*x + A.zy*y + A.zz*z;
        double r2 = rx*rx + ry*ry + rz*rz;
        if (r2 >= a_0*a_0) continue;
        if (r2 != 0.) ff->coeff1 += 1./sqrt(r2);
        double s = r2/(a_0*a_0) - 1.;
        // gamma(r/a_0) = tau(r2/a_0^2 - 1)
        double gam = tau[nu-1];
        for (int k = nu - 2; k >= 0; k--)
          gam = tau[k] + s*gam;
        ff->coeff1 -= gam/a_0;}
  // compute f->coeff2
  double pi = 4.*atan(1.);
  ff->coeff2 = pi * aL * aL / ((4 * nu + 2) * detA );

  // build grid-to-grid stencil for levels L, L + 1
  // ff->khat[L] = d^{L+1} + DFT of khat^L
  Triple gd = *(Triple *)ff->topGridDim;
  // determine range of kernel evaluations khat^L
  int kdmax = 2*ff->nLim + 1;
  Triple sd = gd;
  sd.x = kdmax < gd.x ? kdmax : gd.x;
  sd.y = kdmax < gd.y ? kdmax : gd.y;
  sd.z = kdmax < gd.z ? kdmax : gd.z;
  // calculate level L kappa hat
  double *kh = ff->khat[L];
  msm4g_tic();
  dAL(ff, gd, gd, kh, detA);
  ff->time_stencil[L] = msm4g_toc();

  // build grid-to-grid stencil for levels L-1, ..., 1
  for (int l = L - 1; l > 0; l--){
    gd.x *= 2; gd.y *= 2; gd.z *= 2;
    sd.x = kdmax < gd.x ? kdmax : gd.x;
    sd.y = kdmax < gd.y ? kdmax : gd.y;
    sd.z = kdmax < gd.z ? kdmax : gd.z;
    double *kh = ff->khat[l];
    for (int i = 0; i < sd.x*sd.y*sd.z; i++) kh[i] = 0.;
    //:::o.time = clock();
    msm4g_tic();
    kaphatA(ff, l, gd, sd, kh, as);
    ff->time_stencil[l] = msm4g_toc();
    //:::clock_t end = clock();
    //:::printf("l = %d, elapsed time = %f\n",
    //:::      l, (double)(end - o.time)/CLOCKS_PER_SEC);
  }
}

void FF_delete(FF *ff) {
  free(ff->aCut);
  free(ff->tau);
  free(ff->Q);
  free(ff->J);
  for (int l = 1; l <= ff->maxLevel; l++){
    free(ff->khat[l]);
    for (int alpha = 0; alpha < 3; alpha++)
      free(ff->omegap[l][alpha]);}
  free(ff->khat);
  free(ff->omegap);
  free(ff->dsigma);
  free(ff->nlist);
  for (int alpha = 0; alpha < 3; alpha++)
    free(ff->cL[alpha]);
#ifdef NO_FFT
  free(ff->fftw_in);
#else
  fftw_destroy_plan(ff->forward);
  fftw_destroy_plan(ff->backward);
  fftw_free(ff->fftw_in);
#endif
  free(ff);
}

//helper functions

double invert(Matrix *A) {  // invert A
  // returns |det A|
  double (*a)[3] = (double (*)[3])A;
  // Matrix Ai; double (*ai)[3] = (double (*)[3])&Ai;
  double ai[3][3];
  // compute adjugate
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      ai[i][j] = a[(j+1)%3][(i+1)%3]*a[(j+2)%3][(i+2)%3]
                - a[(j+2)%3][(i+1)%3]*a[(j+1)%3][(i+2)%3];
  // A adj(A) = (det A) I
  double detA = a[0][0]*ai[0][0] + a[0][1]*ai[1][0] + a[0][2]*ai[2][0];
  for (int i = 0; i < 3; i++) 
    for (int j = 0; j < 3; j++)
      a[i][j] = ai[i][j]/detA;
  return fabs(detA);}

void omegap(FF *ff){
  // construct ff->omegap and ff->cL
  int nu = ff->orderAcc;
  // nu positive even integer
  int nuby2 = nu/2;

  // Compute Phi(n), n = 0, 1,..., nuby2-1
  double *Phi = (double *)calloc(nuby2, sizeof(double));
  for(int i = 1; i < nuby2; i++)  Phi[nuby2-i] = ff->Q[i*nu];
  for (int k = 0; k < nu; k++) Phi[0] += ff->Q[(nuby2-1)*nu + k];

  int M[3] = {ff->topGridDim[0], ff->topGridDim[1], ff->topGridDim[2]};
  double pi = 4.*atan(1.);
  ff->omegap = (double *(*)[3])malloc((ff->maxLevel+1)*sizeof(double *[3]));
  int opLim = 2*ff->nLim;
  for (int l = ff->maxLevel; l >= 1; l--){
    for (int alpha = 0; alpha < 3; alpha++){
      int Ma = M[alpha];
      double *c = (double *)malloc(Ma*sizeof(double));
      for (int k = 0; k < Ma; k++){
        c[k] = Phi[0];
        double theta = 2.*pi*(double)k/(double)Ma;
        for (int n = 1; n <= nu/2 - 1; n++)
          c[k] += 2.*Phi[n]*cos((double)n*theta);
        c[k] *= c[k];
        c[k] = 1./c[k];
      }
      int opD = Ma/2 + 1;
      if (Ma/2 > opLim) opD = opLim + 1;
      ff->omegap[l][alpha] = (double *)calloc(opD, sizeof(double));
      for (int m = 0; m < opD; m++){
        double om_m = c[0];
        for (int k = 1; k <= (Ma - 1)/2; k++)
          om_m += 2.*c[k]*cos(2.*pi*(double)(k*m)/(double)Ma);
        if (Ma%2 == 0) om_m += c[Ma/2]*cos(pi*(double)(m));
        ff->omegap[l][alpha][m] = om_m/(double)Ma;}
      if (l == ff->maxLevel) ff->cL[alpha] = c;
      else free(c);
    }
    M[0] *= 2, M[1] *= 2, M[2] *= 2;}
  free(Phi);
}

static void dAL(FF *ff, Triple gd, Triple sd, double kh[], double detA){
  // add d^{L+1}(A) to kh
  Matrix Ai = *(Matrix *)ff-> Ai;
  int nu = ff->orderAcc;
  double *sigmad = ff->dsigma;
  int L = ff->maxLevel;
  double aL = ff->aCut[L];
  double pi = 4.*atan(1.);
  int cx = 2 * ((gd.x+1)/2) + 1;
  int cy = 2 * ((gd.y+1)/2) + 1;
  int cz = (gd.z+1)/2 + 1;
  double ***psi = (double ***)calloc(cx,sizeof(double **));
  for (int i = 0 ; i < cx ;i++) {
    psi[i] = (double **)calloc(cy,sizeof(double *));
    for (int j = 0; j < cy ; j++)
      psi[i][j] = (double *)calloc(cz,sizeof(double));
  }
  // To make negative indexing possible
  for (int i = 0 ; i < cx ;i++)
    psi[i] += (gd.y+1)/2;
  psi += (gd.x+1)/2;
  double twonum1fac = 1.0;
  for (int i = 1 ; i <= 2*nu - 1 ; i += 2)
    twonum1fac  *= i;
  
  // Taylor coefficients
  double tc[4][7] = {
    {105.0,-1890.0,83160.0,-6486480.0,778377600.0,-132324192000.0,30169915776000.0},
    {10395.0,-270270.0, 16216200.0  ,-1654052400.0,251415964800.0,-52797352608000.0,14572069319808000.0},
  {2027025.0,-68918850.0,5237832600.0,-659966907600.0,121433910998400.0,-30358477749600000.0,9836146790870400000.0},
    {654729075.0,-27498621150.0,2529873145800.0,-379480971870000.0,81967889923920000.0,-23770688077936800000.0,8842695964992489600000.0}
  };
  
  psi[0][0][0] = 0.0;
  for (int kx = -(gd.x+ 1)/2; kx <= (gd.x+1)/2; kx++){
    for (int ky = -(gd.y+1)/2; ky <= (gd.y+1)/2; ky++){
      for (int kz = 0; kz <= (gd.z+1)/2; kz++){
        double ax = Ai.xx * kx + Ai.yx * ky + Ai.zx * kz;
        double ay = Ai.xy * kx + Ai.yy * ky + Ai.zy * kz;
        double az = Ai.xz * kx + Ai.yz * ky + Ai.zz * kz;
        double k2 = ax * ax + ay * ay + az * az;
        if (k2 != 0) {
          double k = sqrt(k2);
          double t = pi*k*aL;
          double t2 = t*t;
          double t4 = t2*t2;
          double t6 = t4*t2;
          double t8 = t6*t2;
          double cosfac = 0.0, sinfac = 0.0;
          //double taylor = 0;
          // TODO: remove unnecessary calculations with if-else
          if (nu == 4) {
            sinfac = 15.0 - 6.0 * t2 ;
            cosfac = 15.0 - 1.0 * t2 ;
            //taylor = 1.0/105.0 - t2/1890.0 + t4/83160.0;
          } else if (nu == 6) {
            sinfac = 945.0 - 420.0 * t2 + 15.0 * t4 ;
            cosfac = 945.0 - 105.0 * t2 +  1.0 * t4 ;
            //taylor = 1.0/10395.0 - t2/270270.0 + t4/16216200.0  - t6/1654052400.0  ;
          } else if (nu == 8) {
            sinfac = 135135.0 - 62370.0 * t2 + 3150.0 * t4 - 28.0 * t6 ;
            cosfac = 135135.0 - 17325.0 * t2 +  378.0 * t4 -  1.0 * t6 ;
            //taylor = 1.0/2027025.0 - t2/68918850.0 + t4/5237832600.0 - t6/659966907600.0 + t8/121433910998400.0 ;
          } else if (nu == 10) {
            sinfac = 34459425.0 - 16216200.0 * t2 + 945945.0 * t4 - 13860.0 * t6 + 45.0 * t8 ;
            cosfac = 34459425.0 -  4729725.0 * t2 + 135135.0 * t4 -   990.0 * t6 +  1.0 * t8 ;
            //taylor = 1.0/654729075.0 - t2/27498621150.0 + t4/2529873145800.0 - t6/379480971870000.0 + t8/81967889923920000.0 - t10/23770688077936800000.0 ;
          }
          
          
          double sinvalue = sin(t)/t;
          double cosvalue = cos(t);
          if (t <= ff->tmax*pi/4.0) {
            double taylor = 0.0;
            double tpower = 1.0;
            for (int i = 0 ; i < ff->ntc ; i++) {
              taylor += tpower/tc[nu/2-2][i];
              tpower *= t2;
            }
            psi[kx][ky][kz] = taylor/t2;
          } else {
            psi[kx][ky][kz]  = (sinfac * sinvalue - cosfac * cosvalue) ;
            psi[kx][ky][kz]  = psi[kx][ky][kz] / pow(t,2.0*nu);
          }
          psi[kx][ky][kz] = psi[kx][ky][kz] * (twonum1fac*pi*aL*aL/detA) ;
          
          /*double sum = 0;
          for (int j = nu/2 ; j<=nu -1 ;j++) {
            sum += aL/(k*detA) * pow(-1,j) *
                  (-cos(pi*k*aL) * sigmad[2*j]   / pow(pi*k*aL,2*j+1)
                   +sin(pi*k*aL) * sigmad[2*j+1] / pow(pi*k*aL,2*j+2));
          }
          psi[kx][ky][kz] = sum; */
        }
      }
    }
  }
  
  for (int kx = 0; kx < gd.x; kx++){
    double cLx = ff->cL[0][kx];
    int ckx = kx < (gd.x+1)/2 ? kx : kx - gd.x  ; // kx < Mx/2
    for (int ky = 0; ky < gd.y; ky++){
      double cLxy = cLx*ff->cL[1][ky]; 
      int cky = ky < (gd.y+1)/2.0 ? ky : ky - gd.y ; // ky < My/2
      for (int kz = 0; kz < (gd.z+1)/2 ; kz++){
        double cLxyz = cLxy*ff->cL[2][kz];
        double chi_cL = psi[ckx][cky][kz] * cLxyz;
        kh[(kx*gd.y + ky)*gd.z + kz] += chi_cL*gd.x*gd.y*gd.z;
      }
      for (int kz = (gd.z+1)/2 ; kz < gd.z; kz++){
        double cLxyz = cLxy*ff->cL[2][kz];
        double chi_cL = psi[-ckx][-cky][gd.z-kz] * cLxyz;
        kh[(kx*gd.y + ky)*gd.z + kz] += chi_cL*gd.x*gd.y*gd.z;
      }
    }
  }
  // Back to original offset
  psi -= (gd.x+1)/2;
  for (int i = 0 ; i < cx ;i++)
    psi[i] -= (gd.y+1)/2;
  for (int i=0;i<cx;i++) {
    for (int j = 0 ; j<cy;j++) {
      free(psi[i][j]);
    }
    free(psi[i]);
  }
  free(psi);
}

void grid2grid(FF* restrict ff,const int l,const Triple gd, double* restrict el, const double* restrict ql,
               const Triple sd, const double* restrict kh){
  msm4g_tic();
  int dmax = 2 * ff->nLim + 1;
  int nu = ff->orderAcc;
  int gdxnew = gd.x + max(min(dmax,gd.x),nu);
  int gdynew = gd.y + max(min(dmax,gd.y),nu);
  int gdznew = gd.z + max(min(dmax,gd.z),nu);
  int ofx = 0, ofy = 0, ofz = 0;
  if (nu > min(dmax,gd.x)) ofx = nu/2 - min(dmax,gd.x)/2;
  if (nu > min(dmax,gd.y)) ofy = nu/2 - min(dmax,gd.y)/2;
  if (nu > min(dmax,gd.z)) ofz = nu/2 - min(dmax,gd.z)/2;
 
  for (int ex = 0 ; ex < gd.x ; ex++) {
    for (int kx = 0 ; kx < sd.x ; kx++) {
      for (int ey = 0 ; ey < gd.y ; ey++) {
        for (int ky = 0 ; ky < sd.y ; ky++) {
          for (int ez = 0 ; ez < gd.z ; ez++) {
            for (int kz = 0 ; kz < sd.z ; kz++) {
              int qx = ex + kx ; 
              int qy = ey + ky ;
              int qz = ez + kz ;
              int e_idx = ex * gd.y * gd.z + ey * gd.z + ez; 
              int q_idx = (qx + ofx) * gdynew * gdznew + (qy + ofy) * gdznew + (qz + ofz);
              int k_idx = kx * sd.y * sd.z + ky * sd.z + kz;
              el[e_idx] += kh[k_idx] * ql[q_idx];
            }
          }
        }
      }
    }
  }
  ff->time_grid2grid[l] = msm4g_toc();
}

static void DFT(Triple gd, double dL[], double khatL[]){
  double twopi = 8.*atan(1.);
  for (int kx = 0; kx < gd.x; kx++)
    for (int ky = 0; ky < gd.y; ky++)
      for (int kz = 0; kz < gd.z; kz++){
        double dLk = 0.;
        double cx = twopi*(double)kx/gd.x;
        double cy = twopi*(double)ky/gd.y;
        double cz = twopi*(double)kz/gd.z;
        for (int nx = 0; nx < gd.x; nx++)
          for (int ny = 0.; ny < gd.y; ny++)
            for (int nz = 0.; nz < gd.z; nz++)
              dLk += cos(cx*(double)nx + cy*(double)ny + cz*(double)nz)
                *khatL[(nx*gd.y + ny)*gd.z + nz];
        dL[(kx*gd.y + ky)*gd.z + kz] = dLk;
      }
}

static void FFT(FF *ff, Triple gd, double dL[], double khatL[]){
#ifdef NO_FFT
  ;
#else
  // dL = DFT of khatL
  for (int i = 0; i < gd.x*gd.y*gd.z; i++)
    ff->fftw_in[i] = (fftw_complex)khatL[i];
  fftw_execute(ff->forward);
  for (int i = 0; i < gd.x*gd.y*gd.z; i++)
    dL[i] = creal(ff->fftw_in[i]);
#endif
}

static double kappaA(FF *ff, int l, Vector s, Vector as);  // kappa_l(A s; A)
static void kaphatA(FF *ff, int l, Triple gd, Triple sd, double kh[],
                    Vector as){
  // add kappaHat_l to kh
  //clock_t begin = clock();
  int kdmax = 2*ff->nLim + 1;
  int kdx = kdmax < gd.x ? kdmax : gd.x;
  int kdy = kdmax < gd.y ? kdmax : gd.y;
  int kdz = kdmax < gd.z ? kdmax : gd.z;
  // construct array of kappa values
  double *kap = (double *)malloc(kdx*kdy*kdz*sizeof(double));
  Vector s;
  for (int i = - kdx/2; i <= (kdx - 1)/2; i++){
    double *kapi = kap + (i + kdx/2)*kdy*kdz;
    s.x = (double)i/(double)gd.x;
    for (int j = - kdy/2; j <= (kdy - 1)/2; j++){
      double *kapij = kapi + (j + kdy/2)*kdz;
      s.y = (double)j/(double)gd.y;
      for (int k = - kdz/2; k <= (kdz - 1)/2; k++){
        s.z = (double)k/(double)gd.z;
        kapij[k + kdz/2] = kappaA(ff, l, s, as);}}}
  if (strcmp(o.test, "test_kappaA") == 0)
    for (int i = 0; i < kdx*kdy*kdz; i++) o.kappa[l][i] = kap[i];
  double **op = ff->omegap[l];
  int opLim = 2*ff->nLim;
  //clock_t end = clock();
  //-printf("elapsed time = %f, iterations = %d\n",
  //-(double)(end - o.time)/CLOCKS_PER_SEC, kdx*kdy*kdz);
  //begin = clock();
  //construct kappa hat element by element
    for (int i1 = - sd.x/2; i1 <= (sd.x - 1)/2; i1++){
      double *khi = kh + (i1 + sd.x/2)*sd.y*sd.z;
      for (int j1 = - sd.y/2; j1 <= (sd.y - 1)/2; j1++){
        double *khij = khi + (j1 + sd.y/2)*sd.z;
        for (int k1 = - sd.z/2; k1 <= (sd.z - 1)/2; k1++){
          double khijk = 0.;
          for (int i0 = - kdx/2; i0 <= (kdx - 1)/2; i0++){
            int i = (i0 - i1 + (3*gd.x)/2)%gd.x - gd.x/2;
            double opi = op[0][abs(i)];
            for (int j0 = - kdy/2; j0 <= (kdy - 1)/2; j0++){
              int j = (j0 - j1 + (3*gd.y)/2)%gd.y - gd.y/2;
              double opij = opi*op[1][abs(j)];
              for (int k0 = - kdz/2; k0 <= (kdz - 1)/2; k0++){
                int k = (k0 - k1 + (3*gd.z)/2)%gd.z - gd.z/2;
                double opijk = opij*op[2][abs(k)];
                double kapijk
                  = kap[((i0 + kdx/2)*kdy + j0 + kdy/2)*kdz + k0 + kdz/2];
                khijk += opijk*kapijk;
              }}}
          khij[k1 + sd.z/2] = khijk;}}}
  free(kap);
  //end = clock();
  //-printf("elapsed time = %f, iterations = %d\n",
  //-(double)(end - o.time)/CLOCKS_PER_SEC, sd.x*sd.y*sd.z);
}

static double kappaA(FF *ff, int l, Vector s, Vector as){
  // kappa_l(A s; A)  ::: might separate l < L and l = L
  Matrix A = *(Matrix *)ff->A;
  double *tau = ff->tau;
  int nu = ff->orderAcc;
  double a_l = ff->aCut[l];
  double rootPi = sqrt(4.*atan(1.));

  double pxmin = floor(s.x - a_l*as.x + 1.), pxmax = ceil(s.x + a_l*as.x - 1.);
  double pymin = floor(s.y - a_l*as.y + 1.), pymax = ceil(s.y + a_l*as.y - 1.);
  double pzmin = floor(s.z - a_l*as.z + 1.), pzmax = ceil(s.z + a_l*as.z - 1.);
  double kap_s = 0.;
  for (double px = pxmin; px <= pxmax; px++)
    for (double py = pymin; py <= pymax; py++)
      for (double pz = pzmin; pz <= pzmax; pz++){ // could do just pz >= 0
        double x = s.x - px, y = s.y - py, z = s.z - pz;
        double rx = A.xx*x + A.xy*y + A.xz*z;
        double ry = A.yx*x + A.yy*y + A.yz*z;
        double rz = A.zx*x + A.zy*y + A.zz*z;
        double r = sqrt(rx*rx + ry*ry + rz*rz);
        double rho = r/a_l;
        if (rho >= 1.) continue;
        // one could instead precompute a list of Vectors (p_i, p_j, p_k)
        double gam, gam2;
        if (rho >= 1.) gam = 1./rho;
        else{ // gam = tau(rho^2 - 1)
          double s = rho*rho - 1.;
          gam = tau[nu-1];
          for (int k = nu-2; k >= 0; k--)
            gam = tau[k] + s*gam;}
        gam /= a_l;
        rho *= 2.;
        if (rho >= 1.) gam2 = 1./rho;
        else{
          double s = rho*rho - 1.;
          gam2 = tau[nu-1];
          for (int k = nu-2; k >= 0; k--)
            gam2 = tau[k] + s*gam2;}
        gam2 *= 2./a_l;
        kap_s += gam2 - gam;}
  return kap_s;}

//::: not used
double kappa(FF *ff, int l, double s[3], double edges[3][3]){
  Matrix Ai = *(Matrix *)edges;
  double detA = invert(&Ai);
  *(Matrix *)ff->Ai = Ai;
  Vector as = {sqrt(Ai.xx*Ai.xx + Ai.yx*Ai.yx + Ai.zx*Ai.zx),
               sqrt(Ai.xy*Ai.xy + Ai.yy*Ai.yy + Ai.zy*Ai.zy),
               sqrt(Ai.xz*Ai.xz + Ai.yz*Ai.yz + Ai.zz*Ai.zz)};
  return kappaA(ff, l, *(Vector *)s, as);}

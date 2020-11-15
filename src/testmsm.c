#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "forcefield.h"

struct O o = {""};

int main(int argc, char **argv){
  double edge[3][3];
  FILE *ifile = fopen("../data/Box44.ini", "r");
  for (int i = 0; i < 9; i++){
    fscanf(ifile, "%lf", &edge[0][i]);
  }
  int N; fscanf(ifile, "%d", &N);
  double *q = (double *)malloc(N * sizeof(double));
  double (*r)[3] = (double (*)[3])malloc(N * sizeof(double[3]));
  for (int i = 0; i < N; i++)
    fscanf(ifile, "%lf %lf   %lf     %lf   %*s %*s %*s %*s", \
           q+i, r[i], r[i]+1, r[i]+2);
  fclose(ifile);
  ifile = fopen("../data/Box44.pot", "r");
  double enXact; fscanf(ifile, "%lf", &enXact);
  fclose(ifile);
  ifile = fopen("../data/Box44.acc", "r");
  fscanf(ifile, "%*s") ;
  double (*Fxact)[3] = (double (*)[3])malloc(N * sizeof(double[3]));
  for (int i = 0; i < N; i++) fscanf(ifile, "%lf", &Fxact[i][0]);
  for (int i = 0; i < N; i++) fscanf(ifile, "%lf", &Fxact[i][1]);
  for (int i = 0; i < N; i++) fscanf(ifile, "%lf", &Fxact[i][2]);
  fclose(ifile);
  for (int i = 0; i < N; i++) {
    Fxact[i][0] *= -q[i];
    Fxact[i][1] *= -q[i];
    Fxact[i][2] *= -q[i];
  }
  double Q2 = 0.;
  for (int i = 0; i < N; i++) Q2 += q[i]*q[i];
  double (*F)[3] = (double (*)[3])malloc(N * sizeof(double[3]));


  int      M4[10] = {44   ,  50,   60,  64,  70,  72,  75,  80, 90,108};
  double  a04[10] = {12.74,11.7,10.36,9.92,9.35,9.17,8.93,8.55,7.9,  7};
  int      M6[9]  = {44   , 50, 60,  64,  70,  72,  75,  80,90};
  double  a06[9]  = {10.64,9.6,8.3,7.88,7.34,7.17,6.94,6.59, 6};
  for (int ord = 4; ord <= 10; ord += 2){
    int icount = ord == 4 ? 10 : 9;
    printf("\n order = %d\n", ord);
    printf("%3s %11s %11s %11s %10s %9s %11s\n","M","h","a0","ferror","ttime","etime","error");
    for (int i = 0; i < icount; i++){
      int Mx = (ord == 4 ? M4[i] : M6[i]);
      int M[3] = {Mx, Mx, Mx};
      double a0 = (ord == 4 ? a04[i] : a06[i]);
      double h = edge[0][0]/(double)Mx;
      
      double avg_time_total = 0.0;
      double avg_time_energy = 0.0;
      double min_time_total = 10000;
      double min_time_energy = 10000;
      double avg_ferr = 0.0;
      double avg_perr = 0.0;
      for (int j = 0 ; j < 10 ; j++) {
        FF *ff = FF_new(N,q,edge);
        FF_set_orderAcc(ff, ord);
        FF_set_topGridDim(ff, M);
        FF_set_cutoff(ff, a0);
        clock_t time1 = clock();
        FF_build(ff);
        time1 = clock() - time1;
        double t1 = (double)(time1)/CLOCKS_PER_SEC;
        clock_t time2 = clock();
        double en = FF_energy(ff, F, r, NULL);
        time2 = clock() - time2;
        double t2 = (double)(time2)/CLOCKS_PER_SEC;
        double Ferr = 0.;
        for (int i = 0; i < N; i++){
          double dFx = F[i][0] - Fxact[i][0];
          double dFy = F[i][1] - Fxact[i][1];
          double dFz = F[i][2] - Fxact[i][2];
          double dF = dFx*dFx + dFy*dFy + dFz*dFz;
          Ferr += dF;}
        Ferr = sqrt(Ferr/(double)N);
        double Fref = (Q2/(double)N)/pow(ff->detA/(double)N, 2./3.);
        FF_delete(ff);
        avg_time_total += t1 + t2;
        avg_time_energy += t2 ;
        avg_ferr += Ferr/Fref;
        avg_perr += en - enXact;
        if (t1 + t2 < min_time_total ) min_time_total = t1 + t2;
        if (t2      < min_time_energy) min_time_energy = t2 ;
      }
      printf("%3d %11.8f %11.8f %11.3e %10.3f %9.3f %11g\n",
                    Mx,h, a0, avg_ferr/10, min_time_total, min_time_energy, avg_perr/10);
    }
  }
}

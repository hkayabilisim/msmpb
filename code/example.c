// example.c
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "forcefield.h"
struct O o = {"csCl"}; // struct for optional output
int main(int argc, char** argv){
  int N = 65850;
  // in this example, the length unit is an Angstrom
  double edges[3][3] = {88., 0., 0.,  0., 88., 0.,  0., 0., 88.};
  // in this example, the charge unit is an elementary charge
  double *q = (double *)malloc(N * sizeof(double));
  double (*r)[3] = (double (*)[3])malloc(N * sizeof(double[3]));
  FILE *ifile = fopen("h2o.dat", "r");
  for (int i = 0; i < N; i++)
    fscanf(ifile, "%lf  %lf   %lf     %lf   %*s %*s %*s %*s",
                   q+i, r[i], r[i]+1, r[i]+2);
  fclose(ifile);
  FF *ff = FF_new(N, q);
  // build a Coulomb force field, with 1 Ang margin in neighborlist calculation
  double margin = 0.0;
  FF_set_cutoff(ff,4.0);
  FF_set_maxLevel(ff, 1);
  FF_build(ff,N,r,edges,margin);

  // parameter values can be accessed:
  printf("error estimate = %f\n", FF_get_estErr(ff));
  printf("number of levels = %d\n", FF_get_maxLevel(ff));
  printf("short-range cutoff = %f\n", FF_get_cutoff(ff));
  printf("order of accuracy = %d\n", FF_get_orderAcc(ff));
  printf("using FFT = %s\n", FF_get_FFT(ff) ? "true" : "false");
  int dim[3]; FF_get_topGridDim(ff, dim);
  printf("top level grid is %d x %d x %d\n", dim[0], dim[1], dim[1]);
  printf("number of grid levels = %d\n", FF_get_maxLevel(ff));

  double (*F)[3] = (double (*)[3])malloc(N * sizeof(double[3]));
  printf("energy = %f\n", FF_energy(ff, F, r, NULL));
  FF_delete(ff);
}

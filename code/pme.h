#include <stdbool.h>
#include <time.h>
#include "fftw3.h"

typedef struct Vector {double x, y, z;} Vector;
typedef struct Matrix {double xx, xy, xz, yx, yy, yz, zx, zy, zz;} Matrix;
typedef struct Triple {int x, y, z;} Triple;

struct O {
  char *test;
  clock_t time;
}; // struct for optional output
extern struct O o;

// it is recommended that these values be set with the FF_set methods
typedef struct FF {
  int N;
  double *q;
  double A[3][3], Ai[3][3];  // A inverse transpose
  double detA;
  double cutoff;  // abs cutorr
  int orderAcc;
  int topGridDim[3];
  int *nlist;
  int nlist_len;
  int nlist_interaction_count;
  double tolDir;
  fftw_complex *fftw_in;
  fftw_plan forward, backward;
  double beta;
  double *Q;  // B-spline coefficients
  // first ord/2 pieces of B-splines Q_ord(t)
  // piece_i(t) = q_i0 + q_i1*(t - i) + ... + q_{i,ord-1}*(t - i)^{ord-1}
  double *cL[3]; // c_x^2, c_y^2, c_z^2:
  // coeffs for interpolating reciprocal sum
  double *khat;  // grid2grid stencil; has dimensions topGridDim
  double coeff1, coeff2;  // coeffs of const part of energy
  double time_partcl2partcl;
  double time_nlist;
  double time_grid2grid;
  double time_stencil;
  double time_anterpolation;
  double time_interpolation;
} FF;
FF *FF_new(int N, double *q, double edges[3][3]);
void FF_set_cutoff(FF *ff, double cutoff);
double FF_get_cutoff(FF *ff);
void FF_set_orderAcc(FF *ff, int orderAcc);
void FF_set_topGridDim(FF *ff, int topGridDim[3]);
void FF_set_tolDir(FF *ff, double tolDir);
void FF_build(FF *ff, double (*position)[3]);
double FF_get_cutoff(FF *ff);
int FF_get_orderAcc(FF *ff);
void FF_get_topGridDim(FF *ff, int topGridDim[3]);
double FF_get_tolDir(FF *ff);
double FF_get_errEst(FF *ff, int N, double *charge);
void FF_rebuild(FF *ff, double edges[3][3], double (*position)[3]);
double FF_energy(FF *ff, double (*force)[3], double (*position)[3], double *weight);
  // if weights == NULL, unit weights are assumed; otherwise
  // weights should point to an array of length FF_get_maxLevel(ff) + 1
void FF_delete(FF *ff);

double msm4g_tictocmanager(int push);
void msm4g_tic();
double msm4g_toc();

Vector prod(Matrix m, Vector v);
Vector prodT(Matrix m, Vector v);

#include <stdbool.h>
#include <time.h>
#include "fftw3.h"

typedef struct Vector {double x, y, z;} Vector;
typedef struct Matrix {double xx, xy, xz, yx, yy, yz, zx, zy, zz;} Matrix;
typedef struct Triple {int x, y, z;} Triple;

struct O {
  char *test;
  double *khatLp1;
  double **kappa;
  double e_const;
  double *e; // e[0] = e^0, e[l] = e^{l:}
  bool recX; // e^{L+1} excluded
  clock_t time;
  int level; // 0, 1 (1 to L-1), 2 (L and L+1)
}; // struct for optional output
extern struct O o;

// it is recommended that these values be set with the FF_set methods
typedef struct FF {
  int N;
  double *q;
  double A[3][3], Ai[3][3];  // A inverse transpose
  double detA;
  double relCutoff;
  double errEst;
  double errEstMoore;
  int ntc ;
  int tmax ;
  double estErr; // tol default 1e-3
  int orderAcc;
  int maxLevel;
  int topGridDim[3];
  int *nlist;
  int nlist_len;
  int nlist_interaction_count;
  bool FFT;
  fftw_complex *fftw_in;
  fftw_plan forward, backward;
  double *tau;  // softener coefficients
  double *dsigma; // d^i sigma(1-)
  double **sigmad;
  double *Q;  // B-spline coefficients
  // first ord/2 pieces of B-splines Q_ord(t)
  // piece_i(t) = q_i0 + q_i1*(t - i) + ... + q_{i,ord-1}*(t - i)^{ord-1}
  double *J;  // 2-scale stencil, indexed from -nu/2 thru nu/2
  double *(*omegap)[3];  // quasi-interpolation coefficients
  int nLim; // ceiling(relCutoff - 1)
  double *aCut; // abs cutoffs
  double cutoff;
  double *cL[3]; // c_x^2, c_y^2, c_z^2:
  // coeffs for interpolating reciprocal sum
  double **khat;  // grid2grid stencils
  // khat[L] has dimensions topGridDim
  //khat[l], l < L, has dimensions that are the lesser of  2*nLim + 1
  //  and those of the grid at level l 
  double coeff1, coeff2;  // coeffs of const part of energy
  double time_partcl2partcl;
  double time_nlist;
  double time_grid2grid[10];
  double time_stencil[10];
  double time_anterpolation;
  double time_interpolation;
  double time_restriction[10];
  double time_prolongation[10];
  double time_padding[10];
} FF;
FF *FF_new(int N, double *q);
double FF_get_cutoff(FF *ff);
void FF_set_cutoff(FF *ff, double cutoff);
void FF_set_orderAcc(FF *ff, int orderAcc);
void FF_set_maxLevel(FF *ff, int maxLevel);
void FF_set_topGridDim(FF *ff, int topGridDim[3]);
void FF_set_FFT(FF *ff, bool FFT);
void FF_build(FF *ff,int N, double (*position)[3],double edges[3][3],double margin);
int FF_get_orderAcc(FF *ff);
int FF_get_maxLevel(FF *ff);
void FF_get_topGridDim(FF *ff, int topGridDim[3]);
bool FF_get_FFT(FF *ff);
double FF_get_errEst(FF *ff);
double FF_get_estErr(FF *ff);
void FF_set_estErr(FF *ff,double estErr);
void FF_rebuild(FF *ff, double edges[3][3], double (*position)[3], double margin);
double FF_energy(FF *ff, double (*force)[3], double (*position)[3], double *weight);
  // if weights == NULL, unit weights are assumed; otherwise
  // weights should point to an array of length FF_get_maxLevel(ff) + 1
void FF_delete(FF *ff);
double FF_get_deltaF(FF *ff);

double msm4g_tictocmanager(int push) ;
void msm4g_tic(void);
double msm4g_toc(void);

double invert(Matrix *A);
void omegap(FF *ff);

#define max(a,b) \
({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b); \
  _a > _b ? _a : _b; })

#define min(a,b) \
({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b); \
  _a > _b ? _b : _a; })

Vector prod(Matrix m, Vector v);
Vector prodT(Matrix m, Vector v);

double _FF_get_errEst(FF *ff,int nu);

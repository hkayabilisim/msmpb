#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include "pme.h"

struct O o = {"csCl"}; // struct for optional output

void usage() {
  fprintf(stderr,"Usage: pme dataFile "
          "[--a0 absolute_cutoff] "
          "[--nu accuracyOrder] "
          "[-M grid-spacing] \n"
          "[--tol-dir direct_tolerance]\n"
          "[--bincoor file]\n"
          "[--potfile .pot file]\n"
          "[--accfile .acc file]\n");
  exit(1);
}

int read_binary_coordinate_file(const char *fname, int natoms, double *xyz);

int main(int argc, char **argv){
  if (argc == 1) {
    usage();
  }
  int M[3] = {0, 0, 0};
  double energy;
  double edge[3][3];
  double (*r)[3];
  double (*F)[3];
  double (*acc)[3];
  double *mass;
  double *q;
  char inifile[200],accfile[200],potfile[200],bincoorfile[200];
  bool readbincoor = false;

  double a0 = -1;
  int nu = -1, Mtop=-1;
  double tol_dir = -1.0;
  sprintf(inifile,"%s.ini",argv[1]);
  sprintf(accfile,"%s.acc",argv[1]);
  sprintf(potfile,"%s.pot",argv[1]);
  for (int i = 0 ; i < argc ;i++) {
    if (strcmp(argv[i],"--a0") == 0) {
      a0 = atof(argv[i+1]);
    } else if (strcmp(argv[i],"--nu") == 0) {
      int nu = atoi(argv[i+1]);
    } else if (strcmp(argv[i],"-M") == 0) {
      Mtop = atoi(argv[i+1]);
      M[0] = Mtop; M[1] = Mtop; M[2] = Mtop;
    } else if (strcmp(argv[i],"--tol-dir") == 0) {
      tol_dir = atof(argv[i+1]);
    } else if (strcmp(argv[i],"--potfile") == 0)  {
      sprintf(potfile,"%s",argv[i+1]);
    } else if (strcmp(argv[i],"--accfile") == 0)  {
      sprintf(accfile,"%s",argv[i+1]);
    } else if (strcmp(argv[i],"--bincoor") == 0)  {
      sprintf(bincoorfile,"%s",argv[i+1]);
      readbincoor = true;
    }
  }
  

  FILE *ifile = fopen(inifile, "r");
  
  // read it in by line number and character position
  char line[180];
  
  if (fgets(line, sizeof line, ifile) != NULL) 
  sscanf(line,"%lf%lf%lf%lf%lf%lf%lf%lf%lf",
         edge[0],edge[0]+1,edge[0]+2,
         edge[1],edge[1]+1,edge[1]+2,
         edge[2],edge[2]+1,edge[2]+2);
  int N;
  if (fgets(line, sizeof line, ifile) != NULL)
    sscanf(line,"%d", &N);
  
  q = (double *)calloc(N,sizeof(double));
  r = (double (*)[3])calloc(N,sizeof(double [3]));
  F = (double (*)[3])calloc(N,sizeof(double [3]));
  acc = (double (*)[3])calloc(N,sizeof(double [3]));
  mass = (double *)calloc(N,sizeof(double));
  
  for (int i = 0; i < N; i++) {
    if (fgets(line, sizeof line, ifile) != NULL)
      sscanf(line, "%lf%lf%lf%lf%lf", &q[i], r[i], r[i]+1, r[i]+2,&mass[i]);
  }
  
  fclose(ifile);
 
  if (readbincoor) {
    int errcode = read_binary_coordinate_file(bincoorfile, N, (double *) r);
    if (errcode != 0) {
      fprintf(stderr, "Failed to read binary coordinate file \"%s\"\n",
              bincoorfile);
      exit(1);
    }
  }
  
  FF *ff = FF_new(N,q,edge);
  if (a0 > 0) FF_set_cutoff(ff, a0);
  if (Mtop > 0) FF_set_topGridDim(ff, M);
  if (nu > 0) FF_set_orderAcc(ff, nu);
  if (tol_dir > 0) FF_set_tolDir(ff, tol_dir);

  
  msm4g_tic();
  FF_build(ff, r);
  double time_build = msm4g_toc();
  
  msm4g_tic();
  energy = FF_energy(ff, F, r, NULL);
  double time_energy = msm4g_toc();
  FF_get_topGridDim(ff, M);
  
  double time_manual_sum = 0;
  time_manual_sum += ff->time_partcl2partcl + ff->time_nlist +
                     ff->time_grid2grid + ff->time_stencil +
                     ff->time_anterpolation + ff->time_interpolation;
  printf("{\n");
  printf("\"%s\" : %10.8f,\n","time_partcl2partcl",ff->time_partcl2partcl);
  printf("\"%s\" : %10.8f,\n","time_nlist",ff->time_nlist);
  printf("\"%s\" : %10.8f,\n","time_grid2grid",ff->time_grid2grid);
  printf("\"%s\" : %10.8f,\n","time_stencil",ff->time_stencil);
  printf("\"%s\" : %10.8f,\n","time_anterpolation",ff->time_anterpolation);
  printf("\"%s\" : %10.8f,\n","time_interpolation",ff->time_interpolation);
  printf("\"%s\" : %10.8f,\n","time_build",time_build);
  printf("\"%s\" : %10.8f,\n","time_energy",time_energy);
  printf("\"%s\" : %10.8f,\n","time_longrange",time_energy - ff->time_partcl2partcl);
  printf("\"%s\" : %10.8f,\n","time_other",time_build+time_energy-time_manual_sum);
  printf("\"%s\" : %10.8f,\n","time_total",time_build+time_energy);
  printf("\"%s\" : \"%s\",\n", "data",argv[1]);
  printf("\"%s\" : %d,\n", "nu",FF_get_orderAcc(ff));
  printf("\"%s\" : %d,\n", "nlist_len",ff->nlist_len);
  printf("\"%s\" : %d,\n", "nlist_interaction_count",ff->nlist_interaction_count);
  printf("\"%s\" : %f,\n", "beta",ff->beta);
  printf("\"%s\" : %d,\n", "TopLevelMx",M[0]);
  printf("\"%s\" : %d,\n", "TopLevelMy",M[0]);
  printf("\"%s\" : %d,\n", "TopLevelMz",M[1]);
  printf("\"%s\" : %d,\n", "NumberOfLevels",1);
  printf("\"%s\" : [%5.2f, %5.2f, %5.2f],\n", "Edge row1",
         edge[0][0],edge[0][1],edge[0][2]);
  printf("\"%s\" : [%5.2f, %5.2f, %5.2f],\n", "Edge row2",
         edge[1][0],edge[1][1],edge[1][2]);
  printf("\"%s\" : [%5.2f, %5.2f, %5.2f],\n", "Edge row3",
         edge[2][0],edge[2][1],edge[2][2]);
  printf("\"%s\" : %f,\n", "cutoff",FF_get_cutoff(ff));
  printf("\"%s\" : %d,\n", "nbar",0); // nbar not defined
  printf("\"%s\" : %10.3e,\n", "tol_dir",FF_get_tolDir(ff));
  printf("\"%s\" : %.16e,\n", "utotal",energy);
  
  FILE *afile = fopen(accfile, "r");
  if (afile != NULL) {
    if (fgets(line, sizeof line, afile) != NULL) {}
    for (int i = 0 ; i < N ; i++) {
      if (fgets(line, sizeof line, afile) != NULL) sscanf(line, "%lf", &(acc[i][0])); }
    for (int i = 0 ; i < N ; i++) {
      if (fgets(line, sizeof line, afile) != NULL) sscanf(line, "%lf", &(acc[i][1])); }
    for (int i = 0 ; i < N ; i++) {
      if (fgets(line, sizeof line, afile) != NULL) sscanf(line, "%lf", &(acc[i][2])); }
    
    // Converting acceleration into force
    for (int i = 0 ; i < N ; i++) {
      acc[i][0] *= -q[i];
      acc[i][1] *= -q[i];
      acc[i][2] *= -q[i];
    }

    double Q2 = 0.0;
    for (int i = 0; i < N; i++) Q2 += q[i]*q[i];
    double Fref = (Q2/(double)N)/pow(ff->detA/(double)N, 2./3.);
    double Uref = (Q2/(double)N)/pow(ff->detA/(double)N, 1./3.);
    
    double max_acc = 0.;
    double max_acc_err = 0.;
    double ferror = 0.0;
    for (int i = 0; i < N; i++){
      double acci
      = sqrt(acc[i][0]*acc[i][0] + acc[i][1]*acc[i][1] +	acc[i][2]*acc[i][2]);
      max_acc = fmax(max_acc, acci);
      double errx = acc[i][0] - F[i][0],
      erry	= acc[i][1] - F[i][1],
      errz	= acc[i][2] - F[i][2];
      double err2 = errx*errx + erry*erry + errz*errz;
      double err = sqrt(err2);
      max_acc_err = fmax(max_acc_err, err);
      ferror += err2 ;
    }
    ferror = sqrt(ferror/(double)N); 
    printf("\"%s\" : %25.16f,\n", "detA",ff->detA);
    printf("\"%s\" : %25.16e,\n", "deltaF",ferror);
    printf("\"%s\" : %25.16e,\n", "deltaF/Fref",ferror/Fref);
    printf("\"%s\" : %25.16e,\n", "fref",Fref);
    printf("\"%s\" : %25.16e,\n", "uref",Uref);
    printf("\"%s\" : %25.16e,\n", "Q2",Q2);

    fclose(afile);
  }
  
  double uintra = 0.0;
  for (int i = 0 ; i < N ; i += 3) {
    double *rO = r[i];
    double *rHa = r[i+1];
    double *rHb = r[i+2];
    double distance_OHa = sqrt( (rO[0] - rHa[0]) * (rO[0] - rHa[0])
                          + (rO[1] - rHa[1]) * (rO[1] - rHa[1])
                          + (rO[2] - rHa[2]) * (rO[2] - rHa[2]) );
    double distance_OHb = sqrt( (rO[0] - rHb[0]) * (rO[0] - rHb[0])
                          + (rO[1] - rHb[1]) * (rO[1] - rHb[1])
                          + (rO[2] - rHb[2]) * (rO[2] - rHb[2]) );
    double distance_HH = sqrt( (rHa[0] - rHb[0]) * (rHa[0] - rHb[0])
                           + (rHa[1] - rHb[1]) * (rHa[1] - rHb[1])
                           + (rHa[2] - rHb[2]) * (rHa[2] - rHb[2]) );
    double q_O = q[i];
    double q_Ha = q[i+1];
    double q_Hb = q[i+2];
    double energy_OHa = q_O * q_Ha / distance_OHa;
    double energy_OHb = q_O * q_Hb / distance_OHb;
    double energy_HH = q_Ha * q_Hb / distance_HH;

    double delta = energy_OHa + energy_OHb + energy_HH;
    uintra += delta;
  }

  /*
   * Value namd_energy is set to electrostatic potential
   * needs to be converted from kcal/mol to charge^2/distance.
   */
  printf("\"%s\" : %25.16e,\n", "uintra",uintra);
  printf("\"%s\" : %25.16e,\n", "utotal_kcal_mol",332.0636*(energy-uintra));

  



  FILE *pfile = fopen(potfile, "r");
  
  if (pfile != NULL) {
    double energy_expected=0.0;
    if (fgets(line, sizeof line, pfile) != NULL)
      sscanf(line,"%lf", &energy_expected);
    printf("\"%s\" : %25.16e,\n", "poterror",fabs(energy_expected-energy)/fabs(energy_expected));
    printf("\"%s\" : %25.16f}\n", "potabserror",fabs(energy-energy_expected));

    fclose(pfile);
  }
  
  FILE *fp = fopen("pme.acc","w");
  fprintf(fp,"%d\n",N);
  for (int i=0;i<N;i++) fprintf(fp,"%-25.16f\n",-F[i][0]/q[i]);
  for (int i=0;i<N;i++) fprintf(fp,"%-25.16f\n",-F[i][1]/q[i]);
  for (int i=0;i<N;i++) fprintf(fp,"%-25.16f\n",-F[i][2]/q[i]);
  fclose(fp);
  fp = fopen("pme.pot","w");
  fprintf(fp,"%25.16e\n",energy);
  fclose(fp);
  
  FF_delete(ff);
  free(r);
  free(F);
  free(acc);
  free(q);
  free(mass);

}


/*
 * Read NAMD binary coordinate file.
 * Provide file name (full path to file), number of atoms N, and preallocated
 * array buffer space for 3N coordinates stored in x/y/z form.
 */
int read_binary_coordinate_file(const char *fname, int natoms, double *xyz)
{
  FILE *file = NULL;
  int n = 0;
  if ((file=fopen(fname, "rb")) == NULL) {
    fprintf(stderr, "Unable to open file \"%s\" for reading\n", fname);
    return -1;
  }
  else if (fread(&n, sizeof(int), 1, file) != 1) {
    fprintf(stderr, "Unable to read number of atoms from file \"%s\"\n", fname);
    return -1;
  }
  else if (n != natoms) {
    fprintf(stderr,
            "Expecting %d atoms in file \"%s\" but says it contains %d atoms\n",
            natoms, fname, n);
    return -1;
  }
  else if (fread(xyz, sizeof(double), 3*n, file) != 3*n) {
    fprintf(stderr, "Unable to read %d atom coordinates from file \"%s\"\n",
            n, fname);
    return -1;
  }
  else if (fclose(file) != 0) {
    fprintf(stderr, "Unable to close file \"%s\" after reading\n", fname);
    return -1;
  }
  return 0;
}

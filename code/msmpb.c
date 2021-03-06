#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include "forcefield.h"

struct O o = {"csCl"}; // struct for optional output

void usage() {
  fprintf(stderr,"Usage: msmpb dataFile "
          "[--a0 cutoff] "
          "[--nu accuracyOrder] "
          "[-M grid-spacing] \n"
          "[-L numberOfLevels] "
          "[--perturb]\n"
          "[--repl replx reply replz]\n"
          "[--edge a11 a12 a13 a21 a22 a23 a31 a32 a33]\n"
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
  double edge[3][3], edge_in[3][3];
  bool edge_given = false;
  double (*r)[3];
  double (*F)[3];
  double (*acc)[3];
  double *q;
  double *mass;
  char inifile[200],accfile[200],potfile[200],bincoorfile[200];
  double perturbx = 0.0;
  double perturby = 0.0;
  double perturbz = 0.0;
  bool perturb = false;
  bool readbincoor = false;
  int replx = 1; int reply = 1; int replz = 1;
  double a0 = -1;
  int nu = -1,L=-1,Mtop=-1;
  sprintf(inifile,"%s.ini",argv[1]);
  sprintf(accfile,"%s.acc",argv[1]);
  sprintf(potfile,"%s.pot",argv[1]);
  
  for (int i = 0 ; i < argc ;i++) {
    if (strcmp(argv[i],"--a0") == 0) {
      a0 = atof(argv[i+1]);
    } else if (strcmp(argv[i],"--nu") == 0) {
      nu = atoi(argv[i+1]);
    } else if (strcmp(argv[i],"-M") == 0) {
      Mtop = atoi(argv[i+1]);
      M[0] = Mtop; M[1] = Mtop; M[2] = Mtop;
    } else if (strcmp(argv[i],"-L") == 0) {
      L = atoi(argv[i+1]);
    } else if (strcmp(argv[i],"--perturb") == 0) {
      perturb = true;
    } else if (strcmp(argv[i],"--repl") == 0) {
      replx = atoi(argv[i+1]);
      reply = atoi(argv[i+2]);
      replz = atoi(argv[i+3]);
    } else if (strcmp(argv[i],"--edge") == 0) {
      edge_in[0][0] = atof(argv[i+1]);
      edge_in[0][1] = atof(argv[i+2]);
      edge_in[0][2] = atof(argv[i+3]);
      edge_in[1][0] = atof(argv[i+4]);
      edge_in[1][1] = atof(argv[i+5]);
      edge_in[1][2] = atof(argv[i+6]);
      edge_in[2][0] = atof(argv[i+7]);
      edge_in[2][1] = atof(argv[i+8]);
      edge_in[2][2] = atof(argv[i+9]);
      edge_given = true;
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
  
  if (fgets(line, sizeof line, ifile) != NULL) {
    sscanf(line,"%lf%lf%lf%lf%lf%lf%lf%lf%lf",
      edge[0],edge[0]+1,edge[0]+2,
      edge[1],edge[1]+1,edge[1]+2,
      edge[2],edge[2]+1,edge[2]+2);}
  if (edge_given) memcpy(edge,edge_in, 9*sizeof(double));
  int N;
  if (fgets(line, sizeof line, ifile) != NULL)
    sscanf(line,"%d", &N);
  int Nrep = N*replx*reply*replz;
  q = (double *)calloc(Nrep,sizeof(double));
  r = (double (*)[3])calloc(Nrep,sizeof(double [3]));
  F = (double (*)[3])calloc(Nrep,sizeof(double [3]));
  acc = (double (*)[3])calloc(N,sizeof(double [3]));
  mass = (double *)calloc(Nrep,sizeof(double));

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
  
  // Replication
  int m = N ;
  for (int i = 0; i < replx ; i++) {
    for (int j = 0; j < reply ; j++) {
      for (int k = 0; k < replz ; k++) {
        if (i==0 && j == 0 && k == 0 ) continue;
        for (int n = 0 ; n < N ; n++) {
          r[m][0] = r[n][0] + i * edge[0][0] + j * edge[1][0] + k * edge[2][0];
          r[m][1] = r[n][1] + i * edge[0][1] + j * edge[1][1] + k * edge[2][1];
          r[m][2] = r[n][2] + i * edge[0][2] + j * edge[1][2] + k * edge[2][2];
          mass[m] = mass[n];
          q[m] = q[n];
          m++;
        }
      }
    }
  }
  edge[0][0] *= replx ;
  edge[0][1] *= replx ;
  edge[0][2] *= replx ;
  edge[1][0] *= reply ;
  edge[1][1] *= reply ;
  edge[1][2] *= reply ;
  edge[2][0] *= replz ;
  edge[2][1] *= replz ;
  edge[2][2] *= replz ;
  
  if (perturb) {
    srand(time(0));
    rand(); // skip first random number which seems not random on Mac OS
    perturbx = (2.0 * rand()/(double)RAND_MAX - 1.0 )* edge[0][0];
    perturby = (2.0 * rand()/(double)RAND_MAX - 1.0 )* edge[1][1];
    perturbz = (2.0 * rand()/(double)RAND_MAX - 1.0 )* edge[2][2];
    for (int i = 0 ; i < Nrep ; i++) {
      r[i][0] += perturbx ;
      r[i][1] += perturby ;
      r[i][2] += perturbz ;
    }
  }
  
  o.e = (double *)malloc((L+1)*sizeof(double));
  FF *ff = FF_new(Nrep,q);
  if (L > 0) FF_set_maxLevel(ff, L);
  if (Mtop > 0) FF_set_topGridDim(ff, M);
  if (nu > 0) FF_set_orderAcc(ff, nu);
  if (a0 > 0) FF_set_cutoff(ff,a0);
  msm4g_tic();
  double margin = 0.0;
  FF_build(ff,Nrep,r,edge,margin);
  double time_build = msm4g_toc();
  
  msm4g_tic();
  energy = FF_energy(ff, F, r, NULL);
  double time_energy = msm4g_toc();
  FF_get_topGridDim(ff, M);
  
  double time_manual_sum = 0;
  printf("{\n");
  printf("\"%s\" : %10.8f,\n","time_partcl2partcl",ff->time_partcl2partcl);
  printf("\"%s\" : %10.8f,\n","time_nlist",ff->time_nlist);
  time_manual_sum += ff->time_partcl2partcl + ff->time_nlist;
  //double time_grid2grid_total = 0;
  for (int l = 1 ; l <= ff->maxLevel ; l++) {
    //time_grid2grid_total += ff->time_grid2grid[l];
    printf("\"time_grid2grid_atlevel%d\"        : %10.8f,\n",l,ff->time_grid2grid[l]);
    printf("\"time_padding_atlevel%d\"          : %10.8f,\n",l,ff->time_padding[l]);
    printf("\"time_stencil_atlevel%d\"          : %10.8f,\n",l,ff->time_stencil[l]);
    time_manual_sum += ff->time_grid2grid[l];
    time_manual_sum += ff->time_padding[l];
    time_manual_sum += ff->time_stencil[l];
  }
  //printf("%-30s : %10.8f\n","time_grid2grid_total",time_grid2grid_total);
  //double time_restriction_total = 0;
  for (int l = 2 ; l <= ff->maxLevel ; l++) {
    //time_restriction_total += ff->time_restriction[l];
    printf("\"time_restrictionfrom%dto%d\"       : %10.8f,\n",l-1,l,
           ff->time_restriction[l]);
    time_manual_sum += ff->time_restriction[l];
  }
  //printf("%-30s : %10.8f\n","time_restriction_total",time_restriction_total);
  
  //double time_prolongation_total = 0;
  for (int l = ff->maxLevel-1 ; l >= 1 ; l--) {
    //time_prolongation_total += ff->time_prolongation[l];
    printf("\"time_prolongationfrom%dto%d\"      : %10.8f,\n",l+1,l,
           ff->time_prolongation[l]);
    time_manual_sum += ff->time_prolongation[l];
  }
  //printf("%-30s : %10.8f\n","time_prolongation_total",time_prolongation_total);
  
  printf("\"%s\" : %10.8f,\n","time_anterpolation",ff->time_anterpolation);
  printf("\"%s\" : %10.8f,\n","time_interpolation",ff->time_interpolation);
  time_manual_sum += ff->time_anterpolation + ff->time_interpolation;

  printf("\"%s\" : %10.8f,\n","time_build",time_build);
  printf("\"%s\" : %10.8f,\n","time_energy",time_energy);
  printf("\"%s\" : %10.8f,\n","time_longrange",time_energy - ff->time_partcl2partcl);
  printf("\"%s\" : %10.8f,\n","time_other",time_build+time_energy-time_manual_sum);
  printf("\"%s\" : %10.8f,\n","time_total",time_build+time_energy);
  printf("\"%s\" : \"%s\",\n", "data",argv[1]);
  printf("\"%s\" : %d,\n", "NumberOfLevels",FF_get_maxLevel(ff));
  printf("\"%s\" : %d,\n", "Nrep",Nrep);
  printf("\"%s\" : %d,\n", "N",N);
  printf("\"%s\" : [%d, %d, %d],\n", "repl",replx,reply,replz);
  printf("\"%s\" : %-10.5f,\n","Perturbationx",perturbx);
  printf("\"%s\" : %-10.5f,\n","Perturbationy",perturby);
  printf("\"%s\" : %-10.5f,\n","Perturbationz",perturbz);
  printf("\"%s\" : %d,\n", "nu",FF_get_orderAcc(ff));
  printf("\"%s\" : %d,\n", "ntc",ff->ntc);
  printf("\"%s\" : %d,\n", "nlist_len",ff->nlist_len);
  printf("\"%s\" : %d,\n", "nlist_interaction_count",ff->nlist_interaction_count);
  printf("\"%s\" : %d,\n", "tmax",ff->tmax);
  printf("\"%s\" : %f,\n", "cutoff",FF_get_cutoff(ff));
  printf("\"%s\" : %f,\n", "margin",margin);
  printf("\"%s\" : [%5.2f, %5.2f, %5.2f],\n", "Edge row1",
           edge[0][0],edge[0][1],edge[0][2]);
  printf("\"%s\" : [%5.2f, %5.2f, %5.2f],\n", "Edge row2",
           edge[1][0],edge[1][1],edge[1][2]);
  printf("\"%s\" : [%5.2f, %5.2f, %5.2f],\n", "Edge row3",
           edge[2][0],edge[2][1],edge[2][2]);
  printf("\"%s\" : %d,\n", "TopLevelMx",M[0]);
  printf("\"%s\" : %d,\n", "TopLevelMy",M[1]);
  printf("\"%s\" : %d,\n", "TopLevelMz",M[2]);
  printf("\"%s\" : %.16e,\n", "utotal",energy);
  
  FILE *afile = fopen(accfile, "r");
  if (afile != NULL) {
    if (fgets(line, sizeof line, afile) != NULL) {}
    for (int i = 0 ; i < N ; i++) {
      if (fgets(line, sizeof line, afile)!=NULL) sscanf(line, "%lf", &(acc[i][0])); }
    for (int i = 0 ; i < N ; i++) {
      if (fgets(line, sizeof line, afile)!=NULL) sscanf(line, "%lf", &(acc[i][1])); }
    for (int i = 0 ; i < N ; i++) {
      if (fgets(line, sizeof line, afile)!=NULL) sscanf(line, "%lf", &(acc[i][2])); }
    fclose(afile);
    
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
    double rms_top = 0.0;
    double rms_bottom = 0.0;
    double ferror = 0.0;
    
    int m = 0 ;
    for (int i = 0; i < replx ; i++) {
      for (int j = 0; j < reply ; j++) {
        for (int k = 0; k < replz ; k++) {
          for (int n = 0 ; n < N ; n++) {
            double acci2 = acc[n][0]*acc[n][0] + acc[n][1]*acc[n][1] + acc[n][2]*acc[n][2];
            double acci = sqrt(acci2);
            max_acc = fmax(max_acc, acci);
            double errx = acc[n][0] - F[m][0];
            double erry = acc[n][1] - F[m][1];
            double errz = acc[n][2] - F[m][2];
            double err2 = errx*errx + erry*erry + errz*errz;
            ferror += err2 ;
            double err = sqrt(err2);
            max_acc_err = fmax(max_acc_err, err);
            rms_top += err2/mass[m];
            rms_bottom += acci2/mass[m];
            m++;
          }
        }
      }
    }
    ferror = sqrt(ferror/(double)N);
    double rms = sqrt(rms_top/rms_bottom);
    printf("\"%s\" : %25.16e,\n", "deltaF",ferror);
    printf("\"%s\" : %25.16e,\n", "deltaF/Fref",ferror/Fref);
    printf("\"%s\" : %25.16e,\n", "deltaFest",FF_get_deltaF(ff));
    printf("\"%s\" : %25.16e,\n", "deltaFest/Fref",FF_get_errEst(ff));
    printf("\"%s\" : %25.16e,\n", "fref",Fref);
    printf("\"%s\" : %25.16e,\n", "uref",Uref);
    printf("\"%s\" : %25.16e,\n", "Q2",Q2);
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
    energy_expected *= replx * reply * replz;
    printf("\"%s\" : %25.16e,\n", "poterror",fabs(energy_expected-energy)/fabs(energy_expected));
    printf("\"%s\" : %25.16f\n", "potabserror",fabs(energy-energy_expected));
    fclose(pfile);
  }
  printf("}");
  
  FILE *fp = fopen("msmpb.acc","w");
  fprintf(fp,"%d\n",N);
  for (int i=0;i<Nrep;i++) fprintf(fp,"%-25.16f\n",-F[i][0]/q[i]);
  for (int i=0;i<Nrep;i++) fprintf(fp,"%-25.16f\n",-F[i][1]/q[i]);
  for (int i=0;i<Nrep;i++) fprintf(fp,"%-25.16f\n",-F[i][2]/q[i]);
  fclose(fp);
  fp = fopen("msmpb.pot","w");
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

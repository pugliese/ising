#include "stdlib.h"
#include "time.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "metropolis.h"
#include "lattice.h"
int test_pick(int *lattice,int n, int niter);
int test_correlacion(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int niter, int nsaltos, int ks);
int test_metropolis(int *lattice, int n, float B, float J, float* LUT, float* p_e, int* p_m);

int main(int argc, char **argv) {
  int n = 32;
  int *lattice = malloc(n * n * sizeof(int));
  float prob = 0.5;
  float T = 1;
  float J=0;
  float E;
  int M;
  int niter = 20000;
  float B = 0.1;
  srand(time(NULL));
  float LUT[10];
  for(int i=0;i<5;i++){ // Los posibles valores de los spins de alrededor son 2*i-4 para i=0,..,4 (-4,-2,0,2,4)
    LUT[i] = exp(-(2*(J*(2*i-4)+B))/T);  // Spin positivo
    LUT[i+5] = exp((2*(J*(2*i-4)+B))/T);  // Spin negativo
  }
  M=fill_lattice(lattice, n, prob);
  E=energia_0(lattice,n,J,B);
  printf("E=%f\nM=%d\n", E,M);/*
  for (int i = 0; i < niter; i++) {
    metropolis(lattice, n, B, LUT, EM);
  }
  print_lattice(lattice, n);*/

  //test_pick(lattice,n,niter);
  test_correlacion(lattice, n, B, J, LUT, &E, &M, 100, n*n, n*n);
  /*test_metropolis(lattice, n, B,J, LUT, &E, &M);
  printf("%f   %f\n", E, -n*n*B*tanh(B/T));
  printf("%d\n", M);*/
  free (lattice);
  return 0;
}

int test_pick(int *lattice, int n, int niter){
  int *A = malloc(niter*sizeof(int));
    for (int i=0;i<niter;i++){
      A[i]= pick_site(lattice, n);
    }
    FILE* fp = fopen("Test_pick.txt","a");
    fprintf(fp,"Test de la funciòn pick, con %d iteraciones \n", niter );
    for(int j=0; j< niter-1; j++) {
      fprintf(fp, "%d ,", A[j]);
    }
    fprintf(fp, "%d \n", A[niter-1]);
fclose(fp);
free (A);
}

int test_correlacion(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int niter, int nsaltos, int ks){
  int secs = (unsigned) time(NULL);
  float* CE = (float *) malloc(ks*sizeof(float));
  float* CM = (float *) malloc(ks*sizeof(float));
  for(int k=0;k<ks;k++){
    //printf("Voy a calcular la correlacion\n");
    float *corrs = correlacion(lattice, n, B, J, LUT, p_e, p_m, k+1, niter, nsaltos);
    //printf("Asigne la correlacion\n");
    CE[k] = corrs[0];
    CM[k] = corrs[1];
    free(corrs);
    //printf("Libere la memoria!\n");
  }
  //printf("Calcule todo\n");
  FILE* fp = fopen("test_correlacion.txt","a");
  fprintf(fp, "Test de la funcion de correlacion promediando %d correlaciones, calculadas con %d saltos cada una\nE: ", niter, nsaltos);
  for(int k=0;k<ks-1;k++){
    fprintf(fp, "%f, ", CE[k]);
  }
  fprintf(fp, "%f\n", CE[ks-1]);
  fprintf(fp, "M:");
  for(int k=0;k<ks-1;k++){
    fprintf(fp, "%f, ", CM[k]);
  }
  fprintf(fp, "%f\n", CM[ks-1]);
  secs = (unsigned) time(NULL)-secs;
  int mins = secs/60;
  int horas = mins/60;
  fprintf(fp, "La simulacion tomó %d hs, %d min, %d segs\n\n", horas, mins, secs-60*mins-3600*horas);
  //printf("Escribi todo\n");
  fclose(fp);
  free(CE);
  free(CM);
  return 0;
}

int test_metropolis(int *lattice, int n, float B, float J, float* LUT, float* p_e, int* p_m){
  int res=0;
  print_lattice(lattice, n);
  printf("\n\n");
  for(int i=0;i<10000*n*n;i++){
    res = res+metropolis(lattice,n,B,J,LUT,p_e,p_m);
  }
  print_lattice(lattice, n);
  return res;
}

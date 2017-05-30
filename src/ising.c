#include "stdlib.h"
#include "time.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "metropolis.h"
#include "lattice.h"
int test_pick(int *lattice,int n, int niter);
int test_correlacion(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int ks, int niter, int nsaltos);
int test_metropolis(int *lattice, int n, float B, float J, float* LUT, float* p_e, int* p_m);
int test_vecinos(int n, int idx);
int test_LUT(float* LUT);
float magnet(int *lattice, int n, float p,float T_max, float T_min, int T_pasos,float B,float J, int niter,int k);

int main(int argc, char **argv) {
  int n = 32;
  int *lattice = malloc(n * n * sizeof(int));
  float prob = 0.5;
  float T = 10;
  float J=0;
  float E;
  int M;
  int niter = 2000;
  float B = 0.1;
  //srand(1);
  srand(time(NULL));
  float* LUT =LookUpTable(J,B,T);
  M=fill_lattice(lattice, n, prob);
  E=energia_0(lattice,n,J,B);

  magnet(lattice, n, prob, 3, 1.5, 251  , 1 , J, 1000,20000);

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

int test_correlacion(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int ks,int niter, int nsaltos){
  int secs = (unsigned) time(NULL);
  float* CE = (float *) malloc(ks*sizeof(float));
  float* CM = (float *) malloc(ks*sizeof(float));
  float *corrs = (float *) malloc(2*sizeof(float));
  for(int k=0;k<ks;k++){
    //printf("Voy a calcular la correlacion\n");
    printf("Arranca %d\n", k);
    corrs = correlacion(lattice, n, B, J, LUT, p_e, p_m, k, niter, nsaltos);
    //printf("Asigne la correlacion\n");
    CE[k] = corrs[0];
    CM[k] = corrs[1];
    free(corrs);
    //printf("Libere la memoria!\n");
  }
  //printf("Calcule todo\n");
  FILE *fp = fopen("test_correlacion.txt","a");
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

int test_vecinos(int n, int idx){
  int* lattice = (int *) malloc(n*n*sizeof(int));
  for(int i=0;i<n*n;i++){
    lattice[i] = i;
  }
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      printf("%d ", i*n+j);
    }
    printf("\n");
  }
  printf("El elemento a buscar es el %d\n", idx);
  suma_vecinos(lattice,n,idx);
  free(lattice);
  return 0;
}

int test_LUT(float* LUT){
  for(int i=0;i<5;i++){
    printf("%f  =  ", LUT[i]);
    printf("%f <---------->", 1/LUT[i+5]);
    printf("%f  =  ", 1/LUT[i]);
    printf("%f \n", LUT[i+5]);
  }
  return 0;
}


float magnet(int *lattice, int n, float p,float T_max, float T_min, int T_pasos,float B,float J, int niter,int k){
  float Ti;
  int M = fill_lattice(lattice,n,p);
  float E = energia_0(lattice,n,J,B);
  float* LUT;
  FILE *fp = fopen("Magnetizacion.txt","a");
  float T_step = (T_max-T_min)/(T_pasos-1);
  fprintf(fp, "Magnetizacion para Temperaturas entre %g y %g con %d pasos \n ", T_max,T_min,T_pasos);
  int secs = time(NULL);
  for(Ti=T_min;Ti<=T_max;Ti=Ti+T_step){
    float m_T=0;
    LUT = LookUpTable(J,B,Ti);
    for (int j = 0; j < niter; j++){
      for(int l=0;l<k;l++){
        metropolis(lattice, n,B,J,LUT,&E,&M);
      }
      m_T = m_T+((float)M)/niter;
    }
  fprintf(fp, "%g ,", m_T);
  printf("%g finalizado\n", Ti);
  free(LUT);
  }
  secs = time(NULL)-secs;
  fprintf(fp, "\nEl calculo tomo %d hs, %d min, %d secs\n\n", secs/3600, (secs/60) % 60, secs % 60);
  fclose(fp);
}

#include "stdlib.h"
#include "time.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "metropolis.h"
#include "lattice.h"
#include "ejercicios.h"
int test_pick(int *lattice,int n, int niter);
int test_correlacion(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int ks, int niter, int nsaltos);
int test_metropolis(int *lattice, int n, float B, float J, float* LUT, float* p_e, int* p_m);
int test_vecinos(int n, int idx);
int test_LUT(float* LUT);
float magnet(int *lattice, int n, float p,float T_max, float T_min, int T_pasos,float B,float J, int niter,int k);
int test_pick_2(int *lattice, int n, int niter);

int main(int argc, char **argv) {
  int n = 32;
  int *lattice = malloc(n * n * sizeof(int));
  float prob = 0.5;
  float T =1;
  float J=0.2;
  float E;
  int M;
  int niter = 2000;
  float B = 0;
  //srand(1);
  srand(time(NULL));
  float* LUT =LookUpTable(J,B,T);
  M=fill_lattice(lattice, n, prob);
  E=energia_0(lattice,n,J,B);

  ej_2a(lattice, n,  prob,3,  0.3, 4, 1,  0, 1000,10000 );
//  ej_2c(lattice,n,prob,1,0.1,0.6,100,0,10000,2200);

  //printf("%p\n", (void *) lattice);
  //magnet(lattice, n, prob, 3, 1.5, 251  , 1 , J, 1000,20000);

  ej_2b(lattice, n, 0.5, 0.6, 6, 8*n*n, 500, n*n);
  //int secs = time(NULL);
  //int paso = calc_paso(lattice, n, B, J, LUT, &E, &M, 10*n, 10*n);
  //secs = time(NULL)-secs;
  //  printf("Biseccion: %d en %d min, %d segs\n", paso, secs/60, secs%60);

  int Programa;
  sscanf(argv[1], "%x", &Programa);


  if(Programa == 42){
    float T_max, T_min;
    int T_pasos,k,*lattice2;
    sscanf(argv[2], "%d", &n);
    sscanf(argv[3], "%f", &prob);
    sscanf(argv[4], "%f", &T_min);
    sscanf(argv[5], "%f", &T_max);
    sscanf(argv[6], "%d", &T_pasos);
    sscanf(argv[7], "%f", &B);
    sscanf(argv[8], "%d", &niter);
    sscanf(argv[9], "%d", &k);
    lattice2 = (int *) malloc(n*n*sizeof(int));
    ej_2a(lattice2, n,  prob, T_max,  T_min,  T_pasos, B,  0, niter, k);
    free(lattice2);
  }

  if(Programa == 43){
    float J_max, J_min;
    int m,k,*lattice2;
    sscanf(argv[2], "%d", &n);
    sscanf(argv[3], "%f", &J_min);
    sscanf(argv[4], "%f", &J_max);
    sscanf(argv[5], "%d", &m);
    sscanf(argv[6], "%d", &k);
    sscanf(argv[7], "%d", &niter);
    sscanf(argv[8], "%d", &nsaltos);
    lattice2 = (int *) malloc(n*n*sizeof(int));
    ej_2b(lattice2, n, J_min, J_max, m, k, niter, nsaltos);
    free(lattice2);
  }

  if(Programa == 44){
    float J_max, J_min, T;
    int m,k,*lattice2;
    sscanf(argv[2], "%d", &n);
    sscanf(argv[3], "%f", &T);
    sscanf(argv[4], "%f", &J_min);
    sscanf(argv[5], "%f", &J_max);
    sscanf(argv[6], "%d", &m);
    sscanf(argv[7], "%f", &B);
    sscanf(argv[8], "%d", &niter);
    sscanf(argv[9], "%d", &k);
    lattice2 = (int *) malloc(n*n*sizeof(int));
    ej_2c(lattice, n, 0.5,T, J_min,J_max, m, B,niter,k);
    free(lattice2);
  }


  free(LUT);
  free (lattice);

  return 0;
}

int test_pick(int *lattice, int n, int niter){

  int A;
  FILE* fp = fopen("Test_pick.txt","a");
  fprintf(fp,"Test de la funciòn pick, con %d iteraciones \n", niter );
  for (int i=0;i<niter-1;i++){
    A= pick_site(lattice, n);
    fprintf(fp, "%d ,", A);
  }
  A= pick_site(lattice, n);
  fprintf(fp, "%d\n\n", A);
  fclose(fp);
}

int test_pick_2(int *lattice, int n, int niter){
  int paso=0, dim = n*n, i=0;
  for(int j=0;j<niter;j++){
    paso = pick_site(lattice, n);
    i=i+(paso==dim);
  }
  printf("Iteraciones=%d -> %d\n", niter, i);
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
  float Ti=T_min;
  int M = fill_lattice(lattice,n,p);
  float E = energia_0(lattice,n,J,B);
  float* LUT;
  FILE *fp = fopen("Magnetizacion.txt","a");
  float T_step = (T_max-T_min)/(T_pasos-1);
  fprintf(fp, "Magnetizacion para Temperaturas entre %g y %g con %d pasos con J=%f \n ", T_max,T_min,T_pasos,J);
  int secs = time(NULL);
  float *ener_t= malloc (T_pasos*sizeof(float));
  for(int k=0;k<T_pasos;k++){
    Ti=T_min+k*T_step;
    float m_T=0;
    float e_T=0;
    LUT = LookUpTable(J,B,Ti);
    for (int i=0;i<3000;i++){
      metropolis(lattice, n,B,J,LUT,&E,&M);
    }
    for (int j = 0; j < niter; j++){
      for(int l=0;l<k;l++){
        metropolis(lattice, n,B,J,LUT,&E,&M);
      }
      m_T = m_T+((float)M)/niter;
      e_T = e_T+((float)E)/niter;
    }
  ener_t[k]=e_T;
  fprintf(fp, "%g ,", m_T);
  printf("%g finalizado\n", Ti);
  free(LUT);
  }
  fprintf(fp, "\n Energia para Temperaturas entre %g y %g con %d pasos con J=%f \n ", T_max,T_min,T_pasos,J);
  for(int i=0;i<T_pasos;i++){
    fprintf(fp, "%g ,", ener_t[i]);
    printf("%d finalizado\n", i);
  }
  free(ener_t);
  secs = time(NULL)-secs;
  fprintf(fp, "\nEl calculo tomo %d hs, %d min, %d secs\n\n", secs/3600, (secs/60) % 60, secs % 60);
  fclose(fp);
}

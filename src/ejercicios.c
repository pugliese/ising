#include "ejercicios.h"
#include "metropolis.h"
#include "lattice.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int ej_2b(int *lattice, int n, float J_min, float J_max, int m, int niter, int nsaltos){  // Este J en realidad es J/T
  float *LUT,J, J_step=(J_max-J_min)/(m-1);
  float E;
  int M, paso;
  FILE* fp = fopen("Ejercicio_2b.txt","a");
  fprintf(fp, "Correlacion en funcion de J/T entre %g y %g con %d valores\n", J_min, J_max, m);
  int secs = time(NULL);
  for(int i=0;i<m;i++){
    J = J_min + i*J_step;
    LUT = LookUpTable(J,0,1);
    M = fill_lattice(lattice,n,0.5);
    E = energia_0(lattice,n,J,0);
    /* Termalizo, por ahora avanzo 2000 puntos */
    for(int j=0;j<2000;j++){
      metropolis(lattice, n, 0, J, LUT, &E, &M);
    }
    paso = calc_paso(lattice, n, 0, J, LUT, &E, &M, niter, nsaltos);
    fprintf(fp, "%d, ", paso);
    free(LUT);
  }
  secs = time(NULL)-secs;
  fprintf(fp, "\nEl calculo tomo %d hs, %d min, %d secs\n\n", secs/3600, (secs/60) % 60, secs % 60);
  return 0;
}


int ejercicio_2a(int *lattice, int n, float p,float T_max, float T_min, int T_pasos,float B, float J,int niter,int k){
  float Ti;
  float* LUT;
  FILE *fp = fopen("Ejercicio2a.txt","a");
  float T_step = (T_max-T_min)/(T_pasos-1);
  fprintf(fp, "Energia para Temperaturas entre %g y %g con %d pasos \n ", T_max,T_min,T_pasos);
  int secs = time(NULL);
    for(Ti=T_min;Ti<=T_max;Ti=Ti+T_step){
      fprintf(fp,"Temperatura %f \n", Ti );
      printf("Temperatura %f\n", Ti);
      float* m_e= malloc (k*sizeof(float));
      for (int i=0;i<niter;i++){
        int M = fill_lattice(lattice,n,p);
        float E = energia_0(lattice,n,J,B);
        LUT = LookUpTable(J,B,Ti);
        for (int j = 0; j < k; j++){
          metropolis(lattice, n,B,J,LUT,&E,&M);
          m_e[j] = m_e[j] +((float) M)/niter;
          }
        }
        for (int j=0; j<k ;j++){
          fprintf(fp, "%g ,", m_e[j]);
        }
        fprintf(fp,"\n");
        free(m_e);
        free(LUT);
      }
  secs = time(NULL)-secs;
  fprintf(fp, "\nEl calculo tomo %d hs, %d min, %d secs\n\n", secs/3600, (secs/60) % 60, secs % 60);
  fclose(fp);
  return 0;
}

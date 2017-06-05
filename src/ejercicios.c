#include "ejercicios.h"
#include "metropolis.h"
#include "lattice.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int ej_2b(int *lattice, int n, float J_min, float J_max, int m, int k, int niter, int nsaltos){  // Este J en realidad es J/T
// m es la cantidad de J's y k es la maxima distancia a la que calculamos correlaciones
  float *LUT,J, J_step=(J_max-J_min)/(m-1);
  float E;
  printf("Ejercicio 2.b):\n");
  int M, paso, rango=k;
  //int* paso = malloc(m*sizeof(int));
  FILE* fp = fopen("Ejercicio_2b.txt","a");
  fprintf(fp, "Correlacion en funcion de J/T entre %g y %g con %d valores\n", J_min, J_max, m);
  int secs = time(NULL);
  for(int i=0;i<m;i++){
    J = J_min + i*J_step;
    LUT = LookUpTable(J,0,1);
    M = fill_lattice(lattice,n,0.5);
    E = energia_0(lattice,n,J,0);
    printf("%f Iniciado", J);
    /* Termalizo, por ahora avanzo 2000 puntos */
    for(int j=0;j<2000;j++){
      metropolis(lattice, n, 0, J, LUT, &E, &M);
    }
    printf(" -> Termalizado");
    paso=calc_paso(lattice,n,0,J,LUT,&E,&M,&rango,niter,nsaltos);
    printf(" -> Finalizado: %d\n", paso);
    /*paso[i]= 0;
    while(paso[i]<k || corrs_e[paso[i]]>0.1 || corrs_m[paso[i]]>0.1){ paso[i]++;}
    fprintf(fp, "Correlaciones para J/T=%f\nEnergia: [", J);
    for(int j=0;j<k-1;j++){
      fprintf(fp, "%f, ", corrs_e[j]);
    }
    fprintf(fp, "%f];\nMagnetizacion", corrs_e[k-1]);
    for(int j=0;j<k-1;j++){
      fprintf(fp, "%f, ", corrs_m[j]);
    }
    fprintf(fp, "%f];\n", corrs_m[k-1]);
    fprintf(fp, "El paso optimo es %d con C_e=%f y C_m=%f\n", paso[i], corrs_e[paso[i]], corrs_m[paso[i]]);*/
    fprintf(fp, "%d, ", paso);
    free(LUT);
  }
  /*fprintf(fp, "El vector de pasos optimos es:\n[");
  for(int i=0;i<m-1;i++){
    fprintf(fp, "%d, ", paso[i]);
  }
  fprintf(fp, "%d]\n", paso[m-1]);*/
  secs = time(NULL)-secs;
  fprintf(fp, "\nEl calculo tomo %d hs, %d min, %d secs\n\n", secs/3600, (secs/60) % 60, secs % 60);
  fclose(fp);
  //free(paso);
  return 0;
}

//int ej_2d(int* lattice, int n, )


  float Ti;
int ej_2a(int *lattice, int n, float p,float T_max, float T_min, int T_pasos,float B, float J,int niter,int k){
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
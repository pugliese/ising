#include "ejercicios.h"
#include "metropolis.h"
#include "lattice.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int ej_2a(int *lattice, int n, float p,float T_max, float T_min, int T_pasos,float B, float J,int niter,int k){
  float* LUT;
  float Ti;
  printf("Ejercicio 2a)\n");
  FILE *fp = fopen("Ejercicio2a.txt","a");
  float T_step = (T_max-T_min)/(T_pasos-1);
  fprintf(fp, "Energia para Temperaturas entre %g y %g con %d pasos \n ", T_max,T_min,T_pasos);
  int secs = time(NULL);
    for(Ti=T_min;Ti<=T_max;Ti=Ti+T_step){
      fprintf(fp,"Temperatura %f \n", Ti );
      printf("Temperatura %f\n", Ti);
      float* m_e= malloc (k*sizeof(float));
      for (int l=0;l<k;l++){
        m_e[l]=0;
      }
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
    /* Termalizo, por ahora avanzo 5000 puntos */
    for(int j=0;j<5000;j++){
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
int ej_2c(int *lattice, int n, float p,float T, float J_min,float J_max, int Paso,float B,int niter,int k){
  float J=J_min;
  int M;
  printf("Ejercicio 2c)\n");
  float* LUT, E;
  FILE *fp = fopen("Ejercicio 2_c.txt","a");
  float step = (J_max-J_min)/(Paso-1);
  fprintf(fp, "Simulacion para B/T= %f tomando J/T entre %f y %f, haciendo %d pasos \nMagnetizacion:\n", B,J_min,J_max,Paso);
  int secs = time(NULL);
  float *ener_t= malloc (Paso*sizeof(float));
  for(int m=0;m<Paso;m++){
    J=J_min+m*step;
    M = fill_lattice(lattice,n,p);
    E = energia_0(lattice,n,J,B);
    float m_J=0;
    float e_J=0;
    LUT = LookUpTable(J,B,T);
    for (int i=0;i<5000;i++){       //Termalizacion
      metropolis(lattice, n,B,J,LUT,&E,&M);
    }
    for (int j = 0; j < niter; j++){
      for(int l=0;l<k;l++){
        metropolis(lattice, n,B,J,LUT,&E,&M);
      }
      m_J = m_J+((float)M)/niter;
      e_J = e_J+((float)E)/niter;
    }
  ener_t[m]=e_J;
  fprintf(fp, "%g ,", m_J);
  printf("%g finalizado\n", J);
  free(LUT);
  }
  fprintf(fp, "\nEnergia:\n ");
  for(int i=0;i<Paso;i++){
    fprintf(fp, "%g ,", ener_t[i]);
  }
  free(ener_t);
  secs = time(NULL)-secs;
  fprintf(fp, "\nEl calculo tomo %d hs, %d min, %d secs\n\n", secs/3600, (secs/60) % 60, secs % 60);
  fclose(fp);
  return 0;
}
int ej_2d(int *lattice, int n, int var,float Xmin,float Xmax,int Paso, float X1, float X2,int niter,int k){
// La variable var=1,2,3 indica cual es la variable que será iterada
// El orden numerico de las variables es JBT osea que 1->J, 2->B y 3->T
// Los valores X1 y X2 van, en orden, a las 2 variables restantes
  int M;
  float *LUT, E;
  float J, B, T;
  float Jstep=0, Bstep=0, Tstep=0;
  FILE *fp = fopen("Ejercicio 2_d.txt","a");
  if (var==1){  // Proceso de seleccion de variable; cambia el "header" del txt generado
    J = Xmin;
    Jstep = (Xmax-Xmin)/(Paso-1);
    B = X1;
    T = X2;
    fprintf(fp, "Simulacion para B=%f y T=%f variando J entre %f y %f, haciendo %d pasos\nSe promediaron %d valores para cada magnitud, espaciados %d pasos\nMagnetizacion:\n",X1,X2,Xmin,Xmax,Paso,niter,k);
  }else{
    if (var==2){
      B = Xmin;
      Bstep = (Xmax-Xmin)/(Paso-1);
      J = X1;
      T = X2;
      fprintf(fp, "Simulacion para J=%f y T=%f variando B entre %f y %f, haciendo %d pasos\nSe promediaron %d valores para cada magnitud, espaciados %d pasos\nMagnetizacion:\n",X1,X2,Xmin,Xmax,Paso,niter,k);
    }
    if (var==3){
      T = Xmin;
      Tstep = (Xmax-Xmin)/(Paso-1);
      J = X1;
      B = X2;
      fprintf(fp, "Simulacion para J=%f y B=%f variando T entre %f y %f, haciendo %d pasos\nSe promediaron %d valores para cada magnitud, espaciados %d pasos\nMagnetizacion:\n",X1,X2,Xmin,Xmax,Paso,niter,k);
    }
  }
  printf("Ejercicio 2d)\n");
  int secs = time (NULL);
  float *ener_t= malloc (Paso*sizeof(float));
  for(int m=0;m<Paso;m++){
    J = J+Jstep;  // Dependiendo de var, dos de estas variables tienen
    B = B+Bstep;  // su step igual a cero. Por lo tanto, son constantes
    T = T+Tstep;
    M = fill_lattice(lattice,n,0.5);
    E = energia_0(lattice,n,J,B);
    float m_J=0;
    float e_J=0;
    LUT = LookUpTable(J,B,T);
    for (int i=0;i<5000;i++){                    // Termalizacion
      metropolis(lattice, n,B,J,LUT,&E,&M);
    }
    for (int j = 0; j < niter; j++){
      for(int l=0;l<k;l++){
        metropolis(lattice, n,B,J,LUT,&E,&M);
      }
      m_J = m_J+((float)M)/niter;
      e_J = e_J+((float)E)/niter;
    }
    ener_t[m]=e_J;
    fprintf(fp, "%g ,", m_J);
    printf("%d finalizado\n", m+1);
    free(LUT);
  }
  fprintf(fp, "\nEnergia:\n");
  for(int i=0;i<Paso;i++){
    fprintf(fp, "%g ,", ener_t[i]);
  }
  free(ener_t);
  secs = time(NULL)-secs;
  fprintf(fp, "\nEl calculo tomo %d hs, %d min, %d secs\n\n", secs/3600, (secs/60) % 60, secs % 60);
  fclose(fp);
  return 0;
}


int ej_2e(int *lattice, int n, int var,float Xmin,float Xmax,int Paso, float X1, float X2,int niter,int k){
// La variable var=1,2,3 indica cual es la variable que será iterada
// El orden numerico de las variables es JBT osea que 1->J, 2->B y 3->T
// Los valores X1 y X2 van, en orden, a las 2 variables restantes
  int M;
  float *LUT, *LUT2, E;
  float J, B, T;
  float Jstep=0, Bstep=0, Tstep=0;
  FILE *fp = fopen("Ejercicio 2_d.txt","a");
  if (var==1){  // Proceso de seleccion de variable; cambia el "header" del txt generado
    J = Xmin;
    Jstep = (Xmax-Xmin)/(Paso-1);
    B = X1;
    T = X2;
    fprintf(fp, "Simulacion para B=%f y T=%f variando J entre %f y %f, haciendo %d pasos\nSe promediaron %d valores para cada magnitud, espaciados %d pasos\nMagnetizacion:\n",X1,X2,Xmin,Xmax,Paso,niter,k);
  }else{
    if (var==2){
      B = Xmin;
      Bstep = (Xmax-Xmin)/(Paso-1);
      J = X1;
      T = X2;
      fprintf(fp, "Simulacion para J=%f y T=%f variando B entre %f y %f, haciendo %d pasos\nSe promediaron %d valores para cada magnitud, espaciados %d pasos\nMagnetizacion:\n",X1,X2,Xmin,Xmax,Paso,niter,k);
    }
    if (var==3){
      T = Xmin;
      Tstep = (Xmax-Xmin)/(Paso-1);
      J = X1;
      B = X2;
      fprintf(fp, "Simulacion para J=%f y B=%f variando T entre %f y %f, haciendo %d pasos\nSe promediaron %d valores para cada magnitud, espaciados %d pasos\nMagnetizacion:\n",X1,X2,Xmin,Xmax,Paso,niter,k);
    }
  }
  printf("Ejercicio 2e)\n");
  int secs = time (NULL);
  float *ener_t= malloc (Paso*sizeof(float));
  for(int m=0;m<Paso;m++){
    J = J+Jstep;  // Dependiendo de var, dos de estas variables tienen
    B = B+Bstep;  // su step igual a cero. Por lo tanto, son constantes
    T = T+Tstep;
    M = fill_lattice(lattice,n,0.5);
    E = energia_0_segundos_vecinos(lattice,n,J,B);
    float m_J=0;
    float e_J=0;
    LUT = LookUpTable(J,B,T);
    LUT2 = LookUpTable(-J,0,T);
    for (int i=0;i<5000;i++){                    // Termalizacion
      metropolis_segundos_vecinos(lattice, n,B,J,LUT,LUT2,&E,&M);
    }
    for (int j = 0; j < niter; j++){
      for(int l=0;l<k;l++){
        metropolis_segundos_vecinos(lattice, n,B,J,LUT,LUT2,&E,&M);
      }
      m_J = m_J+((float)M)/niter;
      e_J = e_J+((float)E)/niter;
    }
    ener_t[m]=e_J;
    fprintf(fp, "%g ,", m_J);
    printf("%d finalizado\n", m+1);
    free(LUT);
    free(LUT2);
  }
  fprintf(fp, "\nEnergia:\n");
  for(int i=0;i<Paso;i++){
    fprintf(fp, "%g ,", ener_t[i]);
  }
  free(ener_t);
  secs = time(NULL)-secs;
  fprintf(fp, "\nEl calculo tomo %d hs, %d min, %d secs\n\n", secs/3600, (secs/60) % 60, secs % 60);
  fclose(fp);
  return 0;
}

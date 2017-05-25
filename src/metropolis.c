#include "metropolis.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int metropolis(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m){
  int idx = pick_site(lattice, n);
  int res = flip(lattice, n, B, J, LUT, idx, p_e, p_m);
  return res;
}

int pick_site(int *lattice, int n) {  // (float) (rand ()   /RAND_MAX) n*n
  int res = (int) (((float) rand()/RAND_MAX)*n*n);
  return res;
}


int suma_vecinos(int* lattice, int n,int idx){  // Sumo los spins de los vecinos
  int i = idx/n;    // Fila
  int j = idx%n;    // Columna
  int f_sup = (i+n-1)%n;  // Fila del spin superior (la columna es la misma)
  int f_inf = (i+1)%n;    // Fila del spin inferior (la columna es la misma)
  int c_der = (j+1)%n;    // Columna del spin derecho (la fila es la misma)
  int c_izq = (j+n-1)%n;  // Columna del spin izquierdo (la fila es la misma)
  int res=lattice[f_inf*n+j]+lattice[f_sup*n+j]+lattice[i*n+c_der]+lattice[i*n+c_izq];
  //printf("inf: %d \nsup: %d\nder: %d\nizq: %d\n",lattice[f_inf*n+j],lattice[f_sup*n+j],lattice[i*n+c_der],lattice[i*n+c_izq]);
  return res;
  //printf("%d\n",res );
}

int flip(int *lattice, int n, float B, float J, float* LUT, int idx, float *p_e, int* p_m){
  int res=0;
  int s = suma_vecinos(lattice,n,idx);
/*   Implementacion con Look-Up Table
  En teoria, esto puede reemplazar ToDO lo que esta comentado abajo, hay que cambiar
              float T ---> float* LUT*/

  float proba = LUT[s/2+2+5*(1-lattice[idx])/2];
  // Si el spin es negativo, lo busca en los ultimos 5. Si es positivo, lo busca en los primeros 5 (no suma 5)
// Dentro del subarray, accede a la posicion i=S/2+2 que corresponde a vecinos sumando 2*i-4
// Es más fácil entender esto viendo la definición de la LUT en ising.c
  if(rand()<proba*RAND_MAX){
    res=1;
    lattice[idx]=-lattice[idx];
    *p_e=(*p_e)-2*(J*s+B)*lattice[idx];
    *p_m=(*p_m)+2*lattice[idx];
  }
  return res;

/*
  float deltaE = 2*(s+B)*lattice[idx];
  if (deltaE<=0){
    res=1;
    lattice[idx]=-lattice[idx];
    *p_e=*p_e +deltaE;
    *p_m=*p_m+2*lattice[idx];
  }else{
    int r = rand();
    res = (r<exp(-deltaE/T)*RAND_MAX);    // Si es verdadero (1), acepto
    if(res){
      lattice[idx] = -lattice[idx];
      *p_e=*p_e +deltaE;
      *p_m=*p_m+2*lattice[idx];
    }
  }
  return res;*/
}

float energia_0(int *lattice, int n, float J, float B){
  float res=0; // res va a  ser la energia
  for (int i=0;i<n*n;i++) {
    res=res-lattice[i]*(J*suma_vecinos(lattice,n,i)/2+B);
  }
  return res;
}

float* correlacion(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int k, int niter, int nsaltos){
  int i,j,m,ceros=0;
  float* corrs= (float *) malloc(2*sizeof(float));
  corrs[0]=0;
  corrs[1]=0;
  float MjMjk,EjEjk,Ej,Mj,Ej2,Mj2,Eo,Mo;
  for(i=0;i<niter;i++){
    MjMjk = 0;   // Productos cruzados entre X[j] y X[j+k]
    EjEjk = 0;
    Ej = 0;   // Valor medio de X en j
    Mj = 0;
    Ej2 = 0;  // Valor medio de X² en j
    Mj2 = 0;
    Eo = 0;   // Guarda los valores en cada paso de j
    Mo =0;
    for(j=0;j<nsaltos;j++){
      Eo = *p_e;
      Mo = (float) *p_m;
      Ej = Ej+(*p_e)/nsaltos;
      Mj = Ej+(float) *p_m/nsaltos;
      Ej2 = Ej2+(*p_e)*(*p_e)/nsaltos;
      Mj2 = Mj2+(float) (*p_m)*(*p_m)/nsaltos;
      //printf("Asigne los no XjXjk\n");
      for(m=0;m<k;m++){
        metropolis(lattice,n,B,J,LUT,p_e,p_m); // Avanzo k pasos
      }
      //printf("Avance los %d pasos\n", k);
      EjEjk = EjEjk+Eo*(*p_e)/nsaltos;
      MjMjk = MjMjk+Mo*(float)(*p_m)/nsaltos;
    }
    if(Ej2-Ej*Ej==0 || Mj2-Mj*Mj==0){
      //printf("dio cero\n" );
      ceros++;
    }else{
      corrs[0] = corrs[0]+(EjEjk-Ej*Ej)/(niter*(Ej2-Ej*Ej));
      corrs[1] = corrs[1]+(MjMjk-Mj*Mj)/(niter*(Mj2-Mj*Mj));
    }
  //printf("Puedo calcular una correlacion!!\n");
  }
  printf("Saltaron %d ceros \n", ceros);
  //printf("Puedo calcular un promedio de correlaciones!\n");
  return corrs;
}
/*    FUNCION ALTERNATIVA PARA LA CORRELACION [MODULARIZADA]
float* correlacion(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int k, int niter, int nsaltos){
  float* corrs= (float *) malloc(2*sizeof(float));
  corrs[0]=0;   // Correlacion de la energia
  corrs[1]=0;   // Correlacion de la magnetizacion
  float* guarda;
  for (int i=0;i<niter;i++){
    guarda=correlacion_una_muestra(lattice,n,B,J,LUT,p_e,p_m,k,nsaltos);
    corrs[0] = corrs[0]+guarda[0]/niter;   // Correlacion de la energia
    corrs[1] = corrs[1]+guarda[1]/niter;   // Correlacion de la magnetizacion
  }
  return corrs;
}

float* correlacion_una_muestra(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int k, int nsaltos){
  float* Ef = malloc(nsaltos*sizeof(float));
  float* Ei = malloc(nsaltos*sizeof(float));
  float* Mf = malloc(nsaltos*sizeof(float));  // Defino las magnetizaciones como float
  float* Mi = malloc(nsaltos*sizeof(float));  // para evitar problemas al dividir
  for(int i=0;i<nsaltos;i++){
    Ei[i] = *p_e;
    Mi[i] = (float) *p_m;
    for (int j=0;j<k;j++){
      metropolis(lattice,n,B,J,LUT,p_e,p_m); // Avanzo k pasos
    }
    Ef[i] = *p_e;
    Mf[i] = (float) *p_m;
  }
  float* res = malloc(2*sizeof(float));
  res[0] = coef_corr(Ei, Ef, nsaltos); // Calculo las correlaciones para cada
  res[1] = coef_corr(Mi, Mf, nsaltos); // observable 0->Energia y 1->Magnetizacion
  free(Ef);
  free(Ei);
  free(Mf);
  free(Mi);
  return res;
}

float coef_corr(float Xi, float Xf, int n){
  float cruzado = 0; // Terminos cruzados <XjXjk>
  float X2 = 0;     // Terminos <Xj²>
  float X = 0;     // Terminos <X>
  for (int i=0;i<n;i++){
    cruzado = cruzado + Xi[i]*Xf[i]/n;
    X2 = X2+Xi[i]*Xi[i]/n;
    X = X+Xi[i]/n;
  }
  float numerador = cruzado - X*X;
  float denominador = X2-X*X;
  return (numerador/denominador);
}
*/

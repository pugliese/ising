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
  int dim = n*n;
  int r = rand();
  int res = (int) ((((float) r)/RAND_MAX)*n*n);
  while (res == dim){
    r = rand();
    res = (int) ((((float) r)/RAND_MAX)*n*n);
  }
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
  //if (f_sup<n && f_inf<n && c_der<n && c_izq<n && f_sup>=0 && f_inf>=0 && c_der>=0 && c_izq>=0){}else{printf("Cagamos\n");}
  //printf("inf: %d \nsup: %d\nder: %d\nizq: %d\n",lattice[f_inf*n+j],lattice[f_sup*n+j],lattice[i*n+c_der],lattice[i*n+c_izq]);
  return res;
  //printf("%d\n",res );
}

int flip(int *lattice, int n, float B, float J, float* LUT, int idx, float *p_e, int* p_m){
  int res=0;
  int s = suma_vecinos(lattice,n,idx);
  float proba = LUT[s/2+2+5*(1-lattice[idx])/2];
  // Si el spin es negativo, lo busca en los ultimos 5. Si es positivo, lo busca en los primeros 5 (no suma 5)
// Dentro del subarray, accede a la posicion i=S/2+2 que corresponde a vecinos sumando 2*i-4
// Es más fácil entender esto viendo la definición de la LUT en ising.c
  if(((float)rand())/RAND_MAX<proba){
    res=1;
    lattice[idx]=-lattice[idx];
    *p_e=(*p_e)-2*(J*s+B)*lattice[idx];
    *p_m=(*p_m)+2*lattice[idx];
  }
  return res;
}

float energia_0(int *lattice, int n, float J, float B){
  float res=0; // res va a  ser la energia
  for (int i=0;i<n*n;i++) {
    res=res-lattice[i]*(J*suma_vecinos(lattice,n,i)/2+B);
  }
  return res;
}

float* LookUpTable(float J, float B, float T){
  float* res = malloc(10*sizeof(float));
  for(int i=0;i<5;i++){ // Los posibles valores de los spins de alrededor son 2*i-4 para i=0,..,4 (-4,-2,0,2,4)
    res[i] = exp(-(2*(J*(2*i-4)+B))/T);  // Spin positivo
    res[i+5] = exp((2*(J*(2*i-4)+B))/T);  // Spin negativo
  }
  return res;
}
/*
float* correlacion(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int k, int niter, int nsaltos){
  int i,j,m,ceros=0;
  float* corrs= (float *) malloc(2*sizeof(float));
  corrs[0]=1;
  corrs[1]=1;
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
      Mj = Ej+(float) (*p_m)/nsaltos;
      Ej2 = Ej2+(*p_e)*(*p_e)/nsaltos;
      Mj2 = Mj2+((float) (*p_m)*(*p_m))/nsaltos;
      for(m=0;m<k;m++){
        metropolis(lattice,n,B,J,LUT,p_e,p_m); // Avanzo k pasos
      }
      EjEjk = EjEjk+Eo*(*p_e)/nsaltos;
      MjMjk = MjMjk+Mo*((float)(*p_m))/nsaltos;
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
}*/

int calc_paso(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int niter, int nsaltos){
  float* corrs;
  int k=1;
  float corr_inf_e, corr_sup_e, corr_inf_m, corr_sup_m;
  corr_inf_e = 1; // El inf es siempre k/2 y para k=1 => inf=0
  corr_inf_m = 1; // donde la correlacion es total
  corrs = correlacion(lattice, n, B, J,LUT,p_e,p_m, k, niter, nsaltos);
  corr_sup_e = corrs[0];
  corr_sup_m = corrs[1];
  while (corrs[0]>0.1 || corrs[1]>0.1){  // Equivalente a corr_max = max(corrs[0],corrs[1])>0.1
    k = 2*k;
    free(corrs);
    corrs = correlacion(lattice, n, B, J,LUT,p_e,p_m, k, niter, nsaltos);
    corr_inf_e = corr_sup_e;
    corr_inf_m = corr_sup_m;
    corr_sup_e = corrs[0];
    corr_sup_m = corrs[1];
  }
  free(corrs);
  int inf = k/2, sup = k, med=0;
  printf("Ahora busco entre %d y %d \n", k/2,k);
  while (inf+1<sup){
    // Como ya sali del while anterior, se que max(corr_sup_e, corr_sup_m)>0.1
    // pero max(corr_inf_e, corr_inf_m)>0.1; busco el minimo k con max_corr(k)>0.1
    med = (inf+sup)/2;
    corrs = correlacion(lattice, n, B, J,LUT,p_e,p_m, med, niter, nsaltos);
    if(corrs[0]<0.1 && corrs[1]<0.1){
      corr_sup_m = corrs[1];
      corr_sup_e = corrs[0];
      sup = med;
    }else{
      corr_inf_m = corrs[1];
      corr_inf_e = corrs[0];
      inf = med;
    }
    free(corrs);
  }
  return med;
}

/*    FUNCION ALTERNATIVA PARA LA CORRELACION [MODULARIZADA]*/
float* correlacion(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int k, int niter, int nsaltos){
  float* corrs= (float *) malloc(2*sizeof(float));
  corrs[0]=0;   // Correlacion de la energia
  corrs[1]=0;   // Correlacion de la magnetizacion
  float* guarda;
  for (int i=0;i<niter;i++){
    guarda=correlacion_una_muestra(lattice,n,B,J,LUT,p_e,p_m,k,nsaltos);
    corrs[0] = corrs[0]+guarda[0]/niter;   // Correlacion de la energia
    corrs[1] = corrs[1]+guarda[1]/niter;   // Correlacion de la magnetizacion
    free(guarda);
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

float coef_corr(float* Xi, float* Xf, int n){
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

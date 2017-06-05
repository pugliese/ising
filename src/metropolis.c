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
  /*while (res == dim){
    r = rand();
    res = (int) ((((float) r)/RAND_MAX)*n*n);
  }*/
  return (res%(n*n));
}


int suma_vecinos(int* lattice, int n,int idx){  // Sumo los spins de los vecinos
  int i = idx/n;    // Fila
  int j = idx%n;    // Columna
  int f_sup = (i+n-1)%n;  // Fila del spin superior (la columna es la misma)
  int f_inf = (i+1)%n;    // Fila del spin inferior (la columna es la misma)
  int c_der = (j+1)%n;    // Columna del spin derecho (la fila es la misma)
  int c_izq = (j+n-1)%n;  // Columna del spin izquierdo (la fila es la misma)
  int res=lattice[f_inf*n+j]+lattice[f_sup*n+j]+lattice[i*n+c_der]+lattice[i*n+c_izq];
  return res;
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

int calc_paso(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int* k, int niter, int nsaltos){
  int rango = *k;
  float* corrs_e = (float*) malloc(rango*sizeof(float));
  float* corrs_m = (float*) malloc(rango*sizeof(float));
  graf_corr_2(lattice,n,B,J,LUT,p_e,p_m,corrs_e,corrs_m,rango,niter,nsaltos);
  int paso=0;
  while(paso<rango && (corrs_e[paso]>0.1 || corrs_m[paso]>0.1)){ paso++;}
  while(paso==rango){ // Si me quedo corto, duplico el rango hasta llegar
    printf("No alcanzo (%d) ", rango);
    rango = 2*rango;
    free(corrs_e);
    free(corrs_m);
    corrs_e = (float*) malloc(rango*sizeof(float));
    corrs_m = (float*) malloc(rango*sizeof(float));
    graf_corr_2(lattice,n,0,J,LUT,p_e,p_m,corrs_e,corrs_m,rango,niter,nsaltos);
    paso = 0;
    while(paso<rango && (corrs_e[paso]>0.1 || corrs_m[paso]>0.1)){ paso++;}
  }
  while (paso<rango/2) {rango=rango/2;} // Si me zarpe con el rango, lo achico para el siguiente
  free(corrs_e);
  free(corrs_m);
  *k = rango; // Reemplazo por un mejor rango de busqueda (o no)
  return paso;
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

int graf_corr(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int k, int niter, int nsaltos){
  int i, j;
  int N = nsaltos*k;
  int secs = time(NULL);
  float* E = malloc(N*sizeof(float));
  float* M = malloc(N*sizeof(float));
  float* corrs_E = malloc(k*sizeof(float));
  float* corrs_M = malloc(k*sizeof(float));
  for(i=0;i<k;i++){
    corrs_E[k]=0;
    corrs_M[k]=0;
  }
  for(i=0;i<niter;i++){
    for(j=0;j<N;j++){
      metropolis(lattice,n,B,J,LUT,p_e,p_m);
      E[j]=*p_e;
      M[j]=(float) *p_m;
    }
    for(j=0;j<k;j++){
      corrs_E[j] = corrs_E[j]+coef_corr_k(E, N, j, nsaltos)/niter;
      corrs_M[j] = corrs_M[j]+coef_corr_k(M, N, j, nsaltos)/niter;
    }
    printf("Iteracion %d\n", i+1);
  }
  secs = time(NULL)-secs;
  FILE* fp = fopen("Correlaciones.txt","a");
  fprintf(fp, "Correlaciones hasta %d pasos usando GuardaTodo\n", k);
  for(j=0;j<k-1;j++){
    fprintf(fp, "%f, ", corrs_E[j]);
  }
  fprintf(fp, "%f\n", corrs_E[k-1]);
  for(j=0;j<k-1;j++){
    fprintf(fp, "%f, ", corrs_M[j]);
  }
  fprintf(fp, "%f\nLa simulacion tardo: %dhs, %dmin, %dsegs\n\n", corrs_M[k-1], secs/3600, (secs/60)%60, secs%60);
  fclose(fp);
  free(E);
  free(M);
  free(corrs_M);
  free(corrs_E);
  return 0;
}

float coef_corr_k(float* X, int N, int k, int nsaltos){
// Aca necesitamos N=ancho*nsaltos con N la longitud de X y ancho el
// espacio entre puntos iniciales (obviamente, k<=ancho)
  int i,j;
  int ancho = N/nsaltos;
  float cruzado = 0; // Terminos cruzados <XjXjk>
  float X2 = 0;     // Terminos <Xj²>
  float Xo = 0;     // Terminos <Xj>
  for(i=0;i<N;i=i+ancho){
    cruzado = cruzado+X[i]*X[i+k]/nsaltos;
    X2 = X2+X[i]*X[i]/nsaltos;
    Xo = Xo+X[i]/nsaltos;
  }
  float numerador = cruzado-Xo*Xo;
  float denominador = X2-Xo*Xo;
  return (numerador/denominador);
}


int graf_corr_2(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, float* corrs_e, float* corrs_m, int k, int niter, int nsaltos){
// Los arreglos corrs_x deben tener longitud k y representan corrs_x[k]=Correlacion de la variable X a k+1 pasos
  int i,j;
  float* corrs_muestra_E = malloc(k*sizeof(float));
  float* corrs_muestra_M = malloc(k*sizeof(float));
  for(j=0;j<k;j++){
    corrs_e[j]=0;
    corrs_m[j]=0;
  }
  for(i=0;i<niter;i++){
    correlaciones(lattice,n,B,J,LUT,p_e,p_m,k,nsaltos,corrs_muestra_E,corrs_muestra_M);
    for(j=0;j<k;j++){
      corrs_e[j]=corrs_e[j]+corrs_muestra_E[j]/niter;
      corrs_m[j]=corrs_m[j]+corrs_muestra_M[j]/niter;
    }
  }
  free(corrs_muestra_E);
  free(corrs_muestra_M);
  return 0;
}


int correlaciones(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int k, int nsaltos, float* corrs_e, float* corrs_m){
  int i,j;
  float* cruzados_E = malloc(k*sizeof(float));
  float* cruzados_M = malloc(k*sizeof(float));
  float* Eo = malloc(nsaltos*sizeof(float));
  float* Mo = malloc(nsaltos*sizeof(float));
  for(i=0;i<k;i++){
    cruzados_M[i] = 0;
    cruzados_E[i] = 0;
  }
  for(i=0;i<nsaltos;i++){
    Mo[i] = (float) *p_m; // Parametros en el estado inicial
    Eo[i] = *p_e;
    for(j=0;j<k;j++){ // Voy promediando los terminos cruzados a medida que avanzo con metropolis
      metropolis(lattice,n,B,J,LUT,p_e,p_m);
      cruzados_M[j] = cruzados_M[j]+Mo[i]*(*p_m)/nsaltos;  // Terminos cruzados a j pasos
      cruzados_E[j] = cruzados_E[j]+Eo[i]*(*p_e)/nsaltos;
    }
  }
  float M = 0, E=0, E2=0, M2=0; // Representan los <X> y los <X²>
  for(i=0;i<nsaltos;i++){
    M = M+Mo[i]/nsaltos;
    E = E+Eo[i]/nsaltos;
    E2 = E2+Eo[i]*Eo[i]/nsaltos;
    M2 = M2+Mo[i]*Mo[i]/nsaltos;
  }
  free(Mo);
  free(Eo);
  float denominador_E = E2-E*E;
  float denominador_M = M2-M*M;
  for(i=0;i<k;i++){
    corrs_e[i] = (cruzados_E[i]-E*E)/denominador_E;
    corrs_m[i] = (cruzados_M[i]-M*M)/denominador_M;
  }
  free(cruzados_E);
  free(cruzados_M);
  return 0;
}
/*
int graf_corr_contiguo(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int k, int niter, int nsaltos){
  int i,j;
  int N = k+nsaltos; // Cantidad total de pasos de metropolis
  float* E = malloc(N*sizeof(float));
  float* M = malloc(N*sizeof(float));
  float* corrs_E = malloc(k*sizeof(float));
  float* corrs_M = malloc(k*sizeof(float));
  for(i=0;i<k;i++){
    corrs_E[k]=0;
    corrs_M[k]=0;
  }
  for(i=0;i<niter;i++){
    for(j=0;j<N;j++){
      metropolis(lattice,n,B,J,LUT,p_e,p_m);
      E[j]=*p_e;
      M[j]=(float) *p_m;
    }
    for(j=0;j<k;j++){
      corrs_E[k] = corrs_E[k]+coef_corr_contig(E, N, k+1, nsaltos)/niter;
      corrs_M[k] = corrs_M[k]+coef_corr_contig(M, N, k+1, nsaltos)/niter;
    }
  }
  free(E);
  free(M);
  return 0;
  }

float coef_corr_contig(float* X, int N, int k, int nsaltos){
  int i,j;
  float cruzados = 0;
  float X2 = 0;
  float Xo = 0;
  for(i=0;i<nsaltos;i++){
    cruzados = cruzados + X[i]*X[i+k]/nsaltos;
    X2 = X2 + X[i]*X[i]/nsaltos;
    Xo = Xo + X[i]/nsaltos;
  }
  float res = (cruzados-Xo*Xo)/(X2-Xo*Xo);
  return res;
}
*/

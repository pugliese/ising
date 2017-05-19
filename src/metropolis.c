#include "metropolis.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int metropolis(int *lattice, int n, float B, float T, float *EM) {
  int idx = pick_site(lattice, n);
  int res = flip(lattice, n, B, T, idx, EM);
  return res;
}

int pick_site(int *lattice, int n) {  // (float) (rand ()   /RAND_MAX) n*n
  /*
  int r = rand();
  while(r >= RAND_MAX-(RAND_MAX % (n*n)) ){
    r = rand();
  }
  int res = r % (n*n);
  */
  int res = (int) (((float) rand()/RAND_MAX)*n*n);
  return res;
}

int suma_vecinos(int* lattice, int n, int idx){  // Sumo los spins de los vecinos
  int i = idx/n;    // Fila
  int j = idx%n;    // Columna
  int res=lattice[((i+1)%n)*n+j]+lattice[((i-1+n)%n)*n+j]+lattice[i*n+((j+1)%n)]+lattice[i*n+((j-1+n)%n)];
  return res;
}

int flip(int *lattice, int n, float B, float T, int idx, float *EM) {
  int res;
  int s = suma_vecinos(lattice,n,idx);
  float deltaE = 2*(s+B)*lattice[idx];
  if (deltaE<=0){
    res=1;
    lattice[idx]=-lattice[idx];
    EM[0]=EM[0] +deltaE;
    EM[1]=EM[1]+2*lattice[idx];
  }else{
    int r = rand();
    res = (r<exp(-deltaE/T)*RAND_MAX);    // Si es verdadero (1), acepto
    if(res){
      lattice[idx] = -lattice[idx];
      EM[0]=EM[0] +deltaE;
      EM[1]=EM[1]+2*lattice[idx];
    }
  }
  return res;
}

float energia_0(int *lattice, int n, float B){
  float res=0; // res va a  ser la energia
  for (int i=0;i<n*n;i++) {
    res=res-lattice[i]*(suma_vecinos(lattice,n,i)/2+B);
  }
  return res;
}

float* correlacion(int* lattice, int n, float B, float T, float* EM, int k, int niter, int nsaltos){
  int i,j,m;
  float corrs= (float *) malloc(2*sizeof(float));
  corrs[0]=0;
  corrs[1]=0;
  float MjMjk,EjEjk,Ej,Mj,Ej2,Mj2;
  for(i=0;i<niter;i++){
    MjMjk = 0;   // Productos cruzados entre X[j] y X[j+k]
    EjEjk = 0;
    Ej = 0;   // Valor medio de X en j
    Mj = 0;
    Ej2 = 0;  // Valor medio de X² en j
    Mj2 = 0;
    for(j=0;j<nsaltos;j++){
      Ej = Ej+EM[0]/nsaltos;
      Mj = Ej+EM[1]/nsaltos;
      Ej2 = Ej2+EM[0]*EM[0]/nsaltos;
      Mj2 = Mj2+EM[1]*EM[1]/nsaltos;
      for(m=0;m<k;m++){
        metropolis(lattice,n,B,T,EM); // Avanzo k pasos
      }
      EjEjk = EjEjk+Ej*EM[0]/nsaltos;
      MjMjk = MjMjk+Mj*EM[1]/nsaltos;
    }
    corrs[0] = corrs[0]+(EjEjk-Ej*Ej)/(niter*(Ej2-Ej*Ej));
    corrs[1] = corrs[0]+(MjMjk-Mj*Mj)/(niter*(Mj2-Mj*Mj));
  }
  return corrs;
}

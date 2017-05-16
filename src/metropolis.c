#include "metropolis.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int metropolis(int *lattice, int n, float T, float *EM) {
  int idx = pick_site(lattice, n);
  flip(lattice, n, T, idx, EM);
  return 0;
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

int flip(int *lattice, int n, float T, int idx, float *EM) {
  int s = 0;
  int i, j, res;
  i = idx/n;    // Fila
  j = idx%n;  // Columna
  s=lattice[((i+1)%n)*n+j]+lattice[((i-1+n)%n)*n+j]+lattice[i*n+((j+1)%n)]+lattice[i*n+((j-1+n)%n)];
  float deltaE = 2*lattice[idx]*s;
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

int energia_0(int *lattice, int n, int B){
  int res=0; // res va a  ser la energia
  for (i=0;i<n*n;i++) {
    res=res-lattice[i]*sum_vecinos(lattice,idx);
  }
  res=res-B*fill_lattice(lattice,n,p);
  return res;
}

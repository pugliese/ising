#include "metropolis.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int metropolis(int *lattice, int n, float T) {
int idx;
  idx = pick_site(lattice, n);
  flip(lattice, n, T, idx);
  return 0;
}

int pick_site(int *lattice, int n) {
  int r = rand();
  while(r >= RAND_MAX-(RAND_MAX % (n*n)) ){
    r = rand();
  }
  int res = r % (n*n);
  return res;
}

int flip(int *lattice, int n, float T, int idx) {
  int s = 0;
  int i, j;
  i = idx/n;    // Fila
  j = idx%n;  // Columna
  s=lattice[((i+1)%n)*n+j]+lattice[((i-1+n)%n)*n+j]+lattice[i*n+((j+1)%n)]+lattice[i*n+((j-1+n)%n)];
  float deltaE = 2*lattice[idx]*s;
  int r = rand();
  int res = (r<exp(-deltaE/T)*RAND_MAX);    // Si es verdadero (1), acepto
  if(res){
    lattice[idx] = -lattice[idx];
  }
  return res;
}

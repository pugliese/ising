#include "metropolis.h"

int metropolis(int *lattice, int n, float T) {
  idx = pick_site(lattice, n);
  flip(lattice, n, T, idx);
  return 0;
}

int pick_site(int *lattice, int n) {
  res = rand() % (n*n);
  return res;
}

int flip(int *lattice, int n, float T, int idx) {
  int s = 0;
  if (idx % n != n-1){
    s = s + lattice[idx+1];
  }
  if (idx % n != 0){
    s = s + lattice[idx-1];
  }
  if (idx/n !=0){
    s = s+lattice[idx-n];
  }
  if (idx/n != n-1){
    s = s+lattice[idx+n];
  }
  float deltaE = 2*lattice[idx]*s;
  r = rand();
  if(r<exp(-deltaE/T)){
    lattice[idx] = -lattice[idx];
  }
  return 0;
}

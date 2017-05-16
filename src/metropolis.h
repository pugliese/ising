#ifndef METROPOLIS_H
#define METROPOLIS_H
int metropolis(int *lattice, int n, float T, float *EM);
int pick_site(int *lattice, int n);
int suma_vecinos(int* lattice, int idx);
int flip(int *lattice, int n, float T, int idx, float B float *EM);
int energia_0(int *lattice,int n, int B);
#endif

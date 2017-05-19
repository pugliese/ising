#ifndef METROPOLIS_H
#define METROPOLIS_H
int metropolis(int *lattice, int n, float B, float T, float *EM);
int pick_site(int *lattice, int n);
int suma_vecinos(int* lattice, int n, int idx);
int flip(int *lattice, int n, float B, float T, int idx, float *EM);
float energia_0(int *lattice,int n, float B);
#endif

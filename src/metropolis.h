#ifndef METROPOLIS_H
#define METROPOLIS_H
int metropolis(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m);
int pick_site(int *lattice, int n);
int suma_vecinos(int* lattice,int n,int idx);
int flip(int *lattice, int n, float B, float J, float* LUT, int idx, float *p_e, int* p_m);
float energia_0(int *lattice,int n, float B);
float* correlacion(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int k, int niter, int nsaltos);
#endif

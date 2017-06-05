#ifndef METROPOLIS_H
#define METROPOLIS_H
int metropolis(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m);
int pick_site(int *lattice, int n);
int suma_vecinos(int* lattice,int n,int idx);
int flip(int *lattice, int n, float B, float J, float* LUT, int idx, float *p_e, int* p_m);
float energia_0(int *lattice,int n, float J, float B);
float* correlacion(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int k, int niter, int nsaltos);
int calc_paso(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int* k, int niter, int nsaltos);
float* LookUpTable(float J, float B, float T);

int graf_corr(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int k, int niter, int nsaltos);
float coef_corr_k(float* X, int N, int k, int nsaltos);

int graf_corr_2(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, float* corrs_e, float* corrs_m, int k, int niter, int nsaltos);
int correlaciones(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int k, int nsaltos, float* corrs_e, float* corrs_m);

float* correlacion_una_muestra(int *lattice, int n, float B, float J, float* LUT, float *p_e, int* p_m, int k, int nsaltos);
float coef_corr(float* Xi, float* Xf, int n);
#endif

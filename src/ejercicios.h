#ifndef EJERCICIOS_H
#define EJERCICIOS_H
int ej_2b(int *lattice, int n, float J_min, float J_max, int m, int k, int niter, int nsaltos);
int ej_2a(int *lattice, int n, float p,float T_max, float T_min, int T_pasos,float B, float J,int niter,int k);
int ej_2c(int *lattice, int n, float p,float T_min, float T_max, int Paso,float J,float B,int niter,int k);
int ej_2d(int *lattice, int n, int var,float Xmin,float Xmax,int Paso, float X1, float X2,int niter,int k);
int ej_2e(int *lattice, int n, int var,float Xmin,float Xmax,int Paso, float X1, float X2,int niter,int k);
#endif

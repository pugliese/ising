#include "ejercicios.h"
#include "metropolis.h"
#include "lattice.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


int ej_2b(int *lattice, int n, float J_min, float J_max, int m, int niter, int nsaltos){  // Este J en realidad es J/T
  float *LUT,J, J_step=(J_max-J_min)/(m-1);
  float E;
  int M, paso;
  FILE* fp = fopen("Ejercicio_2b.txt","a");
  fprintf(fp, "Correlacion en funcion de J/T entre %g y %g con %d valores\n", J_min, J_max, m);
  int secs = time(NULL);
  for(int i=0;i<m;i++){
    J = J_min + i*J_step;
    LUT = LookUpTable(J,0,1);
    M = fill_lattice(lattice,n,0.5);
    E = energia_0(lattice,n,J,0);
    /* Termalizo, por ahora avanzo 2000 puntos */
    for(int j=0;j<2000;j++){
      metropolis(lattice, n, 0, J, LUT, &E, &M);
    }
    paso = calc_paso(lattice, n, 0, J, LUT, &E, &M, niter, nsaltos);
    fprintf(fp, "%d, ", paso);
    free(LUT);
  }
  secs = time(NULL)-secs;
  fprintf(fp, "\nEl calculo tomo %d hs, %d min, %d secs\n\n", secs/3600, (secs/60) % 60, secs % 60);
}

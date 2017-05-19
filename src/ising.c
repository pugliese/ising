#include "stdlib.h"
#include "time.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "metropolis.h"
#include "lattice.h"
int test_pick(int *lattice,int n, int niter);

int main(int argc, char **argv) {
  int n = 32;
  int *lattice = malloc(n * n * sizeof(int));
  float prob = 0.5;
  float T = 2.0;
  float EM[2];
  int niter = 20000;
  srand(time(NULL));
  /*
  EM[1]=fill_lattice(lattice, n, prob);
  EM[0]=energia_0(lattice, n,B);
  for (int i = 0; i < niter; i++) {
    metropolis(lattice, n, T, EM);
  }
  print_lattice(lattice, n);
  */
  test_pick(lattice,n,niter);
  free (lattice);
  return 0;
}

int test_pick(int *lattice, int n, int niter){
  int *A = malloc(niter*sizeof(int));
    for (int i=0;i<niter;i++){
      A[i]= pick_site(lattice, n);
    }
    FILE* fp = fopen("Test_pick.txt","a");
    fprintf(fp,"Test de la funciòn pick, con %d iteraciones \n", niter );
    for(int j=0; j< niter-1; j++) {
      fprintf(fp, "%d ,", A[j]);
    }
    fprintf(fp, "%d \n", A[niter-1]);
fclose(fp);
free (A);
}

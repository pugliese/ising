#include "stdlib.h"
#include "time.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "metropolis.h"
#include "lattice.h"

int main(int argc, char **argv) {
  int n = 32;
  int *lattice = malloc(n * n * sizeof(int));
  float prob = 0.5;
  float T = 2.0;
  float EM[2];
  int niter = 2000;
  srand(time(NULL));
  EM[1]=fill_lattice(lattice, n, prob);
  EM[0]=energia_0(lattice, n,B);
  for (int i = 0; i < niter; i++) {
    metropolis(lattice, n, T, EM);
  }
  print_lattice(lattice, n);
  free (lattice);
  return 0;
}

#include "lattice.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


int fill_lattice(int *lattice, int n, float p) {
  int i,res=0, r;
  for (i=0;i<n*n;i++){
    r= rand();
    lattice[i]=((r<p*RAND_MAX)*2-1);    // Cuando sea verdaddes el <  ==> 1*2-1 ==> 1  // Cuando sea falso el < ==> 0*2-1 ==> -1
    res=res+lattice[i];
  }
  return res;                           //Va a ser la magnetizaciòn osea la cantidad e 1 y -1 de la red
}

int print_lattice(int *lattice, int n) {
  int i,j;
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      if (lattice[i*n+j]==1){
        printf("+ ");
      }else{
        printf("- ");
      }
    }
    printf("\n");
  }
  return 0;
}

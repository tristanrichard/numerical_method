#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "BC.h"

void init_Gaussian(double* U, double* x,double N,double sigma){
    for (int i = 0; i < N; i++){
        U[i]= exp(-pow(x[i]/sigma,2));
    }
}
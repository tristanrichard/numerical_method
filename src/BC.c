#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "BC.h"

void init_Gaussian(double* U, double* x,double N,double sigma,double L,double c){
    for (int i = 0; i < N; i++){
        double a=3.0/5.0;
        U[i]= ((1-a*cos(2*M_PI*(x[i]/L)))/c)*exp(-pow(x[i]/sigma,2));
    }
}
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "BC.h"

void init_Gaussian(double* U, double* x,double N,double sigma,double L,double c,int Uni){
    if(Uni){
        for (int i = 0; i < N; i++){
            double a=3.0/5.0;
            double eta= x[i] - ((a*L)/(2*M_PI))*sin(2*M_PI*(x[i]/L));
            U[i]= ((1-a*cos(2*M_PI*(x[i]/L)))/c)*exp(-pow(eta/sigma,2));
        }
    }else{
        for (int i = 0; i < N; i++){
            U[i]= exp(-pow(x[i]/sigma,2));
        }
    }
}

void init_Wavepacket(double* U, double* x,double N,double sigma,double L,double c,int Uni){
    if(Uni){
        for (int i = 0; i < N; i++){
            double a=3.0/5.0;
            double kp=(2.0*M_PI*12.0)/L;
            double eta= x[i] - ((a*L)/(2*M_PI))*sin(2*M_PI*(x[i]/L));
            U[i]= ((1-a*cos(2*M_PI*(x[i]/L)))/c)*cos(eta*kp)*exp(-pow(eta/sigma,2));
        }
    }else{
        for (int i = 0; i < N; i++){
            double kp=(2.0*M_PI*12.0)/L;
            U[i]= cos(x[i]*kp)*exp(-pow(x[i]/sigma,2));
        }
    }
}
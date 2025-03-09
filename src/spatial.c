#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spatial.h"

int circ(int i, double N){
    return fmod(i+N,N);
}
void (*choose_solver(int scheme))(double* ,double* ,double ,double ,double ,double){
    if(scheme==1){
        return E2;
    }
    else if(scheme==2){
        return E4;
    }else{
        return NULL;
    }
    //else if(scheme==3){
        //return I4;
    //}else{
        //return ED;
    //}
}
void E2(double* U,double* dU,double dt,double N,double c,double h){
    for(int i=0;i<N;i++){
        dU[i] = -dt*c*((U[circ(i+1,N)]-U[circ(i-1,N)])/(2*h));
    }
}

void E4(double* U,double* dU,double dt,double N,double c,double h){
    for(int i=0;i<N;i++){
        double coef1 = (U[circ(i+1,N)]-U[circ(i-1,N)])*(4/3);
        double coef2 = (U[circ(i+2,N)]-U[circ(i-2,N)])*(1/6);
        dU[i] = -dt*c*((coef1-coef2)/(2*h));
    }
}


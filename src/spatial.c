#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spatial.h"
#include "thomas.h"

int circ(int i, double N){
    return fmod(i+N,N);
}
void (*choose_solver(int scheme))(double* ,double* ,double ,double ,double ,double){
    if(scheme==1){
        return E2;
    }
    else if(scheme==2){
        return E4;
    }
    else if(scheme==3){
        return I4;
    }else if(scheme==4){
        return ED;
    }else{
        return NULL;
    }
}
void E2(double* U,double* dU,double dt,double N,double c,double h){
    for(int i=0;i<N;i++){
        dU[i] = -dt*c*((U[circ(i+1,N)]-U[circ(i-1,N)])/(2*h));
    }
}

void E4(double* U,double* dU,double dt,double N,double c,double h){
    for(int i=0;i<N;i++){
        double coef1 = (U[circ(i+1,N)]-U[circ(i-1,N)])*(4.0/3.0);
        double coef2 = (U[circ(i+2,N)]-U[circ(i-2,N)])*(1.0/6.0);
        dU[i] = -dt*c*((coef1-coef2)/(2*h));
    }
}

void I4(double* U,double* dU,double dt,double N,double c,double h){
    double *q = (double *)malloc(N * sizeof(double));
    double a = 3.0/2.0;
    double alpha = 1.0/2.0;
    for(int i=0;i<N;i++){
        q[i] =-dt*c* a *((U[circ(i+1,N)]- U[circ(i-1,N)])/(2*h));
    }
    solve_period_3diag(N,1,alpha/2,alpha/2,dU,q);
    free(q);
}

void ED(double* U,double* dU,double dt,double N,double c,double h){
    for(int i=0;i<N;i++){
        dU[i] = -dt*c*((U[circ(i-2,N)] - (6*U[circ(i-1,N)] + (3*U[circ(i,N)] + (2*U[circ(i+1,N)]))))/(6*h));
    }
}


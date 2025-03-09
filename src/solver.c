#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solver.h"
#include "spatial.h"



void RK4(double* U, double* t, double dt, double c,double N,double h,int scheme) {
    double alpha[4] = {0.0, 0.5, 0.5, 1.0};
    double gamma[4] = {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0};
    double *Us = (double *)malloc(N * sizeof(double));
    double *Ul = (double *)malloc(N * sizeof(double));
    double *delta_u = (double *)malloc(N * sizeof(double));
    //double ts;
    //double tl;//not use here since no dependence on t in f

    void (*solver)(double* ,double* ,double ,double ,double ,double) = choose_solver(scheme);

    for(int i=0;i<N;i++){
        Us[i]=U[i];
        delta_u[i]=0;
    }
    //ts=*t;

    for (int i = 0; i < 4; i++) {
        for(int j=0;j<N;j++){
            Ul[j] = Us[j] + alpha[i]*delta_u[j];
        }
        //tl = ts + alpha[i]*dt;

        solver(U,delta_u,dt,N,c,h);///coder ca doit modifier delta_u

        for(int j=0;j<N;j++){
            U[j] = U[j]+ gamma[i]*delta_u[j];
        }
        *t = *t + gamma[i]*dt;
    }

    free(Us);
    free(Ul);
    free(delta_u);

}
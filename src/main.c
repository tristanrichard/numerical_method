#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "BC.h"
#include "solver.h"

int main(int argc, char *argv[]){
/**
 * main of the code, run all
 * 
 * @param1 scheme= E2/E4/I4/ED 
 * @param2 N = Number of grid points
 * @param3 CFL = CFL number 
 * @param4 nt = number of time steps
 * 
 * @return 0 of OK
 */
    if(argc == 5){
        printf("\n ######################## \n");
        printf("Numerical simulation has started\n");
    }else{
        printf("EXECUTION STOP: Not enough argument");
        return 1;
    }

    //scheme attribution //
    int scheme_number = 0;
    if(strcmp(argv[1],"E2")){
        scheme_number =1;
        printf("Scheme used is E2\n");
    }else if(strcmp(argv[1],"E4")){
        scheme_number= 2;
        printf("Scheme used is E4\n");
    }else{
        printf("Scheme not supported, simulation stopped\n");
        return 1;
    }

    //grid configuration //
    double N = atof(argv[2]);
    double L =1;
    double sigma = L/16;
    double h = L/N;
    double c = 1;
    double *x = (double*)malloc(N * sizeof(double)); //xi
    double CFL = atof(argv[3]);
    double dt = (CFL*h)/c;
    double nt = atof(argv[4]);
    if (x == NULL) {
        printf("Memory allocation failed\n");
        return 1;
    }
    for (int i = 0; i < N; i++) {
        x[i] = -L/2 + i * h;  // Calcul des positions xi
    }
    printf("####### Grid configuration#######\n");
    printf("Number of grid points : %f\n",N);
    printf("Length of the domain: %f\n",L);
    printf("grid spacing: %f\n",h);
    printf("CFL: %f\n",CFL);
    printf("dt: %f\n",dt);
    printf("nt: %f\n",nt);
    

    // Init_condition //
    double *U = (double *)malloc(N * sizeof(double));
    if (U == NULL) {
        printf("Memory allocation failed\n");
        return 1;
    }
    init_Gaussian(U,x,N,sigma);
    printf("####Initialisation sucessfull#####\n");
    printf("U points:");
    for (int i = 0; i < N; i++){
        printf("%f ",U[i]);
    }
    printf("\n");

    // RK4 iteration //
    printf("#### Start of RK4 iteration #####\n");
    double t = 0;
    for(int i=0;i<=nt;i++){
        RK4(U,&t,dt,c,N,h,scheme_number);
        printf("#### End of iteration: %i ####",i);
        printf("t= %f\n",t);
        printf("U points:");
        for (int j = 0; j < N; j++){
            printf("%f ",U[j]);
        }
        printf("\n");
    }




    free(U);
    free(x);

    


}
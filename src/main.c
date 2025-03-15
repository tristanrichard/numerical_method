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
 * @param5 no_uniform si =1
 * @param6 wavepacket si =1
 * 
 * @return 0 of OK
 */
    if(argc == 7){
        printf("\n ######################## \n");
        printf("Numerical simulation has started\n");
    }else{
        printf("EXECUTION STOP: Not enough argument");
        return 1;
    }

    //scheme attribution //
    int scheme_number = 0;
    if(strcmp(argv[1],"E2")==0){
        scheme_number =1;
        printf("Scheme used is E2\n");
    }else if(strcmp(argv[1],"E4")==0){
        scheme_number= 2;
        printf("Scheme used is E4\n");
    }else if(strcmp(argv[1],"I4")==0){
        scheme_number= 3;
        printf("Scheme used is I4\n");
    }else if(strcmp(argv[1],"ED")==0){
        scheme_number= 4;
        printf("Scheme used is ED\n");
    }else{
        printf("Scheme not supported, simulation stopped\n");
        return 1;
    }

    //grid configuration //
    double N = atof(argv[2]);
    double L =1;
    double sigma = L/16;
    double h = L/N;
    double c = 1.0;
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

    //allocation of the convection speed//
    double *b = (double*)malloc(N * sizeof(double));
    if(atof(argv[5])==1){//non uniform case
        for(int i=0;i<N;i++){
            double a = 3.0/5.0;
            b[i]=c/(1-a*cos(2*M_PI*(x[i]/L)));
        }
    }else{
        for(int i=0;i<N;i++){
            b[i]=c;
        }
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

    if(atof(argv[6])==1){
        init_Wavepacket(U,x,N,sigma,L,c,atoi(argv[5]));
    }else{
        init_Gaussian(U,x,N,sigma,L,c,atoi(argv[5]));
    }
    printf("####Initialisation sucessfull#####\n");
    printf("U points:");
    for (int i = 0; i < N; i++){
        printf("%f ",U[i]);
    }
    printf("\n");

    // Writing data //
    char filename[100];  
    sprintf(filename, "data/results_%s_N%.0f_nonuniform%.0f_wave%.0f.txt", argv[1],N,atof(argv[5]),atof(argv[6]));
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file for writing\n");
        return 1;
    }
    fprintf(file, "Grid configuration:\n");
    fprintf(file, "N: %f, L: %f, h: %f, CFL: %f, dt: %f, nt: %f, c: %f, sigma: %f\n", N, L, h, CFL, dt, nt,c,sigma);
    fprintf(file, "Initial condition (U points):\n");
    for (int i = 0; i < N; i++) {
        fprintf(file, "%f ", U[i]);
    }
    fprintf(file, "\n");

    // RK4 iteration //
    printf("#### Start of RK4 iteration #####\n");
    double t = 0;
    for(int i=0;i<=nt;i++){
        RK4(U,&t,dt,b,N,h,scheme_number);
        fprintf(file, "#### End of iteration: %i ####\n", i);
        fprintf(file, "t= %f\n", t);
        fprintf(file, "U points:\n");
        for (int j = 0; j < N; j++) {
            fprintf(file, "%f ", U[j]);
        }
        fprintf(file, "\n");
      
    }
    printf("#### RK4 iterations ended successfully ####");

    fclose(file);


    free(U);
    free(x);
    free(b);

    return 0;

}
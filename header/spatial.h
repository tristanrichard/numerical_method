#ifndef spatial_H
#define spatial_H

int circ(int i, double N);
void (*choose_solver(int scheme))(double* ,double* ,double ,double ,double* ,double);

void E2(double* U,double* dU,double dt,double N,double* c,double h);

void E4(double* U,double* dU,double dt,double N,double* c,double h);

void I4(double* U,double* dU,double dt,double N,double* c,double h);

void ED(double* U,double* dU,double dt,double N,double* c,double h);

#endif
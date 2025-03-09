#ifndef thomas_h
#define thomas_h

void solve_Ac_thomas(const int n, const double a, const double b, const double c,
    double* x,double* q);

void solve_period_3diag(const int n, const double a, const double b, const double c,
    double* x,double* q);

#endif /* thomas_h */
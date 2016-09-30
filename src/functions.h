#ifndef FUNCTIONS_H
#define FUNCTIONS_H

void jacobi(double ** A, double ** R, int n, double epsilon);
double max_offdiag(double ** A, int * k, int * l, int n);
void rotate(double ** A, double ** R, int k, int l, int n);

#endif // FUNCTIONS_H

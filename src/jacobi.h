#ifndef JACOBI_H
#define JACOBI_H
#include <armadillo>
using namespace arma;

void jacobi_method(mat **A, mat **R, int n);

void maxoffdiag(mat **A, int *k, int *l, int n);

#endif // JACOBI_H

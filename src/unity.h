#ifndef UNITY_H
#define UNITY_H

void unity_eig(int *test_counter, double epsilon);
void unity_max(int *test_counter);
double unity_ortho(double ** R, int n, double epsilon);
void unity_init_ortho(double ** A, double ** R, double * V, double *rho,
                double h_marked, int n, double epsilon, int *test_counter);
#endif // UNITY_H

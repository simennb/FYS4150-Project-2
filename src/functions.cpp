#include "functions.h"
#include "unity.h"
#include <cmath>
#include <iostream>
#include <armadillo>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;
using namespace arma;

void jacobi(double ** A, double ** R, int n, double epsilon){
    /*
    Jacobis method
    */

    /* Setting up the eigenvector matrix */
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            if (i == j)
            {
                R[i][j] = 1.0;
            }
            else
            {
                R[i][j] = 0.0;
            }
        }
    }

    /* Initializing some variables */
    int k, l;
    double max_iterations = (double) n * (double) n * (double) n;
    int iterations = 0;
    double maxoffdiag = max_offdiag(A, &k, &l, n);

    /* Rotates the matrix A as long as it has any non-diagonal non-zero elements*/
    while ( fabs(maxoffdiag) > epsilon && (double) iterations < max_iterations){
        rotate(A, R, k, l, n);
        maxoffdiag = max_offdiag(A, &k, &l, n);
        iterations++;
    }

    cout << "Number of iterations " << iterations<< endl;
    ofstream ofile;
    ofile.open("../benchmarks/iterations_n"+to_string(n)+".dat");
    ofile<<setiosflags(ios::showpoint | ios::uppercase);
    ofile<<"Iterations= "<<setw(5)<<iterations;
    ofile.close();

    return;
}



double max_offdiag(double ** A, int * k, int * l, int n){
    /*
    Choosing the maximum nondiagonal element
    setting abs(a[k][l]) = max(abs[i][i+1])
    Assuming that A is symetric round the diagonal
    and therefore the second loop goes only through
    j+1.
    */

    double max = 0.0;
    for (int i = 0; i < n; i++){
        for (int j = i+1; j < n; j++) {
            if(fabs(A[i][j]) > max){
                max = fabs(A[i][j]);
                *l = i;
                *k = j;
            }
        }
    }
    return max;
}

void rotate(double ** A, double ** R, int k, int l, int n){
    /*
    Rotation function. Used to find cos and sin.
    */
    double s, c; //sin and cos functions
    if( A[k][l] != 0.0) {
        double t, tau;
        tau = (A[l][l] - A[k][k])/(2*A[k][l]);
        if (tau > 0) {
            t = -tau + sqrt(1.0 + tau*tau); //want to keep t as small as possible
        }
        else {
            t = -tau - sqrt(1.0 + tau*tau); //if tau is negative, we must subtract t and visa versa
        }
        c = 1.0/sqrt(1 + t*t);
        s = t*c;
    }
    else{
        c = 1.0;
        s = 0.0;
    }

    double a_kk, a_ll, a_il, a_ik, r_ik, r_il;
    a_kk = A[k][k];
    a_ll = A[l][l];

    //finidng diagonal elements
    A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
    A[l][l] = s*s*a_kk + 2.0*s*c*A[k][l] + c*c*a_ll;
    A[k][l] = 0;
    A[l][k] = 0;

    //finding non-diagonal elements
    for (int i = 0; i <n; i++){
        if(i != k && i != l) {
            a_ik = A[i][k];
            a_il = A[i][l];
            A[i][k] = c*a_ik - s*a_il;
            A[k][i] = A[i][k];
            A[i][l] = c*a_il + s*a_ik;
            A[l][i] = A[i][l];
        }
        //new eigenvectors

        r_ik = R[i][k];
        r_il = R[i][l];
        R[i][k] = c*r_ik - s*r_il;
        R[i][l] = c*r_il + s*r_ik;
    }
    return;
}

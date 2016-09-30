#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cmath>
#include <armadillo>
#include <time.h>
#include "functions.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[])
{
    if (argc < 4){
        cout<<"Usage: "<<argv[0]<<" dim"<<" rho_max"<<" task"<<endl;
        cout<<" dim : dimensionality of matrix A"<<endl;
        cout<<" rho_max : boundary for rho"<<endl;
        cout<<" task : noninteract, interact, or unit_test"<<endl;
        exit(1);
    }
    /* Initializing some variables */
    int n = atoi(argv[1]);
    double rho_max = atof(argv[2]);
    double epsilon = 1.0e-8;
    int iterations = 0;

    /* For timing segments of code */
    clock_t start, finish;
    double time_jacobi, time_arma;

    /* Initial conditions */
    double rho_0 = 0.0;
    double h = (rho_max-rho_0)/((double)n);

    /* Create matrices and arrays */
    double **A = new double*[n];  /* matrix */
    double **R = new double*[n];  /* eigenvector matrix */
    for (int i=0; i<n; i++)
    {
        A[i] = new double[n];
        R[i] = new double[n];
    }

    double *V = new double[n];    /* potential array */

    /* Setting up the eigenvector matrix */
    for (int i=0; i<n; i++)
    {
        for (int j=0; i<n; i++)
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

    /* Switches determining which potential to use */
    /////////////////////////////////////////////
    /////////////////////////////////////////////

    // atm just using the simple potential, fix later
    for (int i=1; i<n+1; i++)
    {
        V[i-1] = pow((rho_0 + i*h),2);
    }

    /* Filling matrix A with values */
    double h_marked = 1.0/pow(h,2);

    for (int i=0; i<n-1; i++)
    {
        A[i][i] = 2.0*h_marked + V[i];  /* diagonal elements */
        A[i][i+1] = -h_marked;          /* non diagonal elements */
        A[i+1][i] = -h_marked;
    }
    A[n-1][n-1] = 2.0*h_marked + V[n-1];

    /* Running Jacobi's algorithm */
    start = clock();
    /////////////////////////////////////////////
    jacobi(A, R, n, epsilon);
    /////////////////////////////////////////////
    finish = clock();

    time_jacobi = ((finish-start)/((double)CLOCKS_PER_SEC));

    /* Armadillo solving */
    start = clock();
    // maybe make function to create matrix to use or something
    /////////////////////////////////////////////
    /////////////////////////////////////////////
    finish = clock();

    time_arma = ((finish-start)/((double)CLOCKS_PER_SEC));


    /////////////////////////////////////////////
    /* Interacting case */
    /////////////////////////////////////////////
    if (strcmp(argv[3], "interact") == 0)
    {
        /*for (int i=0; i<n; i++)
        {
            V[i]
        }*/
    }

    for (int i = 0; i < n; i++){
        cout << A[i][i] << endl;
    }
    /* Freeing up memory */
    for (int i=0; i<n; i++)
    {
        delete[] A[i];
        delete[] R[i];
    }
    delete[] A;
    delete[] R;
    delete[] V;

}

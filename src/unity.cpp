#include "unity.h"
#include "functions.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

void unity_eig(int * test_counter, double epsilon){
    /*
    Test to see if our Jacobi method returns the right eigenvalues.
    Should give us eigenvalues 1 and 3
     */
    /*Initializing some variables*/
    int n = 2;
    int eig1 = 1; //known eigenvalues
    int eig2 = 3;
    double **AU = new double*[n];
    double **RU = new double*[n];

    /*For error reports*/
    ofstream unit_ofile;
    unit_ofile.open("../benchmarks/unittest.txt", ios_base::app);

    /*Creating 2x2 matrices*/
    for (int i=0; i<n; i++){
        AU[i] = new double[n];
        RU[i] = new double[n];
    }
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            if (i == j)
            {
                AU[i][j] = 2;
                RU[i][j] = 1;
            }
            else
            {
                AU[i][j] = 1;
                RU[i][j] = 0;
            }
        }
    }

    jacobi(AU, RU, n, epsilon);

    /* Eigenvalue test*/
    if (fabs(AU[0][0] - eig1) < epsilon && fabs(AU[1][1] - eig2) < epsilon)
    {
        *test_counter += 1; //adds to test_counter if test succeded.
    }
    else
    {
        cout << "Eigenvalues test failed"<< endl;
        unit_ofile << "Eigenvalues test failed"<< endl;
    }
    return;
}

void unity_max(int * test_counter)
{   /* Test to see if our function max_offdiag returns
    the maximum non-diagonal value*/
    /* Initializing some variables */
    int n = 5;
    double max_value = 3;
    int k,l;
    double **AU = new double*[n];

    ofstream unit_ofile;
    unit_ofile.open("../benchmarks/unittest.txt",ios_base::app);

    /* Create matrices */
    for (int i=0; i<n; i++)
    {
        AU[i] = new double[n];
    }
    for (int i=0; i<n; i++)
    {
       for (int j=0; j<n;j++)
       {
           if (i == j)
           {
               AU[i][j] = 1;
           }
           else if(i==2 && j == 3)
           {
               //preserve symmetry
               AU[i][j] = max_value;
               AU[j][i] = max_value;
           }
           else
           {
               AU[i][j] = 0;
           }
       }
    }

    double maxoffdiag = max_offdiag(AU, &k, &l, n);

    /* Test to see if the known max value equals the programs
    max value.*/
    if (maxoffdiag == max_value)
    {
        *test_counter += 1;
    }
    else
    {
        cout << "Maximum non-diagonal test failed" << endl;
        unit_ofile << "Maximum non-diagonal test failed" << endl;
    }
    return;
}

double unity_ortho(double ** R, int n, double epsilon)
{
    /* Test to see if a matrix is orthonogal
    Returns the fraction of failed cases. This should
    be equal to zero.*/

    /* Initializing some variables */
    double r = 0;
    double r1 = 0;
    double r2 = 0;
    double fail_counter = 0;
    double counter = 0;
    double failfrac;

    /* Going through the matrix and finding the scalar product
    of an eigenvector with its transposed counterpart.
    R[i][j] gives us the eigenvector column j with eigenvalue at
    point i.*/
    for (int j = 0; j < n; j++)
    {
        for (int k = j; k < n; k++)
        {
            r = 0;
            for (int i = 0; i < n; i++)
            {
                r += R[i][j]*R[i][k]; //scalar product
            }

        if (j == k && fabs(1 - r) < epsilon)
        {
            r1 += r;
        }
        else if (j != k && fabs(0 - r) < epsilon)
        {
            r2 += r;
        }        else
        {
            fail_counter += 1;
        }
        counter += 1;
        }
    }
    failfrac = fabs(0 - fail_counter/counter);
    return failfrac;
}

void unity_init_ortho(double ** A, double ** R, double * V, double *rho,
                double h_marked, int n, double epsilon,
                int * test_counter)
{
    /* Initializes the matrices and vectors needed for
    the non-interacting case to see if we still have orthogonality
    after the rotation.*/
    /////////////////////////////////////////////
    ///         Non-interacting case          ///
    /////////////////////////////////////////////
    for (int i=0; i<n; i++)
    {
        V[i] = pow(rho[i], 2);
    }

    /* Filling matrix A with values */
    for (int i=0; i<n-1; i++)
    {
        A[i][i] = 2.0*h_marked + V[i];  /* diagonal elements */
        A[i][i+1] = -h_marked;          /* non diagonal elements */
        A[i+1][i] = -h_marked;
    }
    A[n-1][n-1] = 2.0*h_marked + V[n-1];

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
    int counter = 0;
    int tot_fail_frac = 0;
    double fail_frac = 0;
    double max_iterations = (double) n * (double) n * (double) n;
    int iterations = 0;
    double maxoffdiag = max_offdiag(A, &k, &l, n);

    ofstream unit_ofile;
    unit_ofile.open("../benchmarks/unittest.txt", ios_base::app);

    /* Jacobi's while-loop with counter*/
    while ( fabs(maxoffdiag) > epsilon && (double) iterations < max_iterations){
        rotate(A, R, k, l, n);
        maxoffdiag = max_offdiag(A, &k, &l, n);
        if (counter == max_iterations/(10*n)) //test of orthogonality
        {
            fail_frac = unity_ortho(R, n, epsilon);
            if (fail_frac < epsilon)
            {
                tot_fail_frac += fail_frac;
            }
            counter = 0;
        }
        counter++;
        iterations++;
    }
    if (fabs(tot_fail_frac) < epsilon)
    {
        *test_counter += 1;
    }
    else
    {
        cout << "Orthogonality test failed" << endl;
        unit_ofile << "Orthogonality test failed" << endl;
    }
    cout << "Number of iterations " << iterations<< endl;
    return;
}

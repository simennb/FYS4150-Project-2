#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cmath>
#include <armadillo>
#include <time.h>
#include "functions.h"
#include "unity.h"

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
    if (strcmp(argv[3], "unit_tests") == 0)
    {
        /* If running unit tests, making sure that no matrices/arrays initialized */
        //unit_testing(); or something when functions made
        exit(1);
    }
    /* Initializing some variables */
    int n = atoi(argv[1]);
    double rho_max = atof(argv[2]);
    double rho_0 = 0.0;
    double h = (rho_max-rho_0)/((double)n); /* step length */
    double h_marked = 1.0/pow(h,2);  /* 1/h^2 is fairly useful */
    double epsilon = 1.0e-8;

    /* For timing segments of code */
    clock_t start, finish;
    double time_jacobi, time_arma;
    ofstream time_ofile;
    time_ofile.open("../benchmarks/timing.txt",ios_base::app);

    /* Create matrices and arrays */
    double **A = new double*[n];  /* matrix */
    double **R = new double*[n];  /* eigenvector matrix */
    for (int i=0; i<n; i++)
    {
        A[i] = new double[n];
        R[i] = new double[n];
    }

    double *V = new double[n];    /* potential array */
    double *rho = new double[n];  /* rho array */

    /* Setting up the eigenvector matrix */
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            if (i == j)
            {
                R[i][i] = 1.0;
            }
            else
            {
                R[i][j] = 0.0;
            }
        }
    }

    /* Initializing rho-array */
    for (int i=1; i<n+1; i++)
    {
        rho[i-1] = rho_0 + i*h; /* unsure about indexes */
    }

    /////////////////////////////////////////////
    ///         Non-interacting case          ///
    /////////////////////////////////////////////
    if (strcmp(argv[3], "non-interact") == 0)
    {
        /* Setting the potential in the non-interacting case */
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

        /* Running Jacobi's algorithm */
        start = clock();
        jacobi(A, R, n, epsilon);
        finish = clock();

        time_jacobi = ((finish-start)/((double)CLOCKS_PER_SEC));

        /* Armadillo solving */
        mat A_arma  = zeros<mat>(n,n);
        vec V_arma  = zeros<vec>(n);
        mat eigvecs = zeros<mat>(n,n);
        vec eigvals = zeros<vec>(n);

        /* Setting values of A_arma to A and V_arma to V */
        for (int i=0; i<n-1; i++)
        {
            A_arma(i,i) = A[i][i];
            A_arma(i,i+1) = A[i][i+1];
            A_arma(i+1,i) = A[i+1][i];
            V_arma(i) = V[i];
        }
        A_arma(n-1,n-1) = A[n-1][n-1];
        V_arma(n-1) = V[n-1];

        start = clock();
        eig_sym(eigvals, eigvecs, A_arma);
        finish = clock();

        time_arma = ((finish-start)/((double)CLOCKS_PER_SEC));
        //for (int i=0; i<3; i++)
        //{
        //   cout<<eigvals[i]<<endl;
        //}

        /* Writing to file */
        ofstream ofile;
        ofile.open("../benchmarks/eigenvalues_noint_n"+to_string(n)+".dat");
        ofile<<setiosflags(ios::showpoint | ios::uppercase);
        ofile<<"N= "<<setw(5)<<n;
        ofile<<" rho_min= "<<setw(3)<<rho_0<<" rho_max= "<<setw(3)<<rho_max<<endl;
        ofile<<"  eigvals jacobi"<<"  eigvals armadillo"<<endl;
        for (int i=0; i<5; i++) /* write first 5 eigenvalues to file */
        {
            ofile<<setw(15)<<setprecision(8)<<A[i][i]<<setw(15)<<setprecision(8)<<eigvals[i]<<endl;
        }
        ofile.close();

        /* Writing execution time to file*/
        time_ofile<<"N= "<<setw(5)<<n;
        time_ofile<<" rho_min= "<<setw(3)<<rho_0<<" rho_max= "<<setw(3)<<rho_max<<endl;
        time_ofile<<"t_jacobi = "<<setw(10)<<setprecision(8)<<time_jacobi;
        time_ofile<<" t_arma = "<<setw(10)<<setprecision(8)<<time_arma<<"\n\n";
        time_ofile.close();
    }

    /////////////////////////////////////////////
    ///           Interacting case            ///
    /////////////////////////////////////////////
    if (strcmp(argv[3], "interact") == 0)
    {
        double omega_r[4] = {0.01, 0.5, 1.0, 5.0}; // small enough to not worry about memory allocation

        for (int i=0; i<4; i++) /* going through the different omega_r values */
        {
            /* Setting the potential in the interacting case */
            for (int j=0; j<n; j++)
            {
                V[j] = pow((rho[j]*omega_r[i]), 2) - 1.0/rho[j];
            }

            /* Filling matrix A with values */
            for (int i=0; i<n-1; i++)
            {
                A[i][i] = 2.0*h_marked + V[i];  /* diagonal elements */
                A[i][i+1] = -h_marked;          /* non diagonal elements */
                A[i+1][i] = -h_marked;
            }
            A[n-1][n-1] = 2.0*h_marked + V[n-1];

            /* Running Jacobi's algorithm */
            jacobi(A, R, n, epsilon);

            /* Writing results to file */
            //ofstream ofile;


        }
    }
    /////////////////////////////////////////////
    ///              Unity tests              ///
    /////////////////////////////////////////////
    unity_eig(epsilon);
    unity_max();

    /////////////////////////////////////////////
    ///            Writing to file            ///
    /////////////////////////////////////////////

    ofstream outfile;



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

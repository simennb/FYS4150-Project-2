#include "unity.h"
#include "functions.h"
#include <iostream>
#include <cmath>
using namespace std;

void unity_eig(double epsilon){
    /*
    Test should give eigenvalues 1 and 3
     */
    int n = 2;
    int eig1 = 1;
    int eig2 = 3;

    double **AU = new double*[n];
    double **RU = new double*[n];

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

    if (fabs(AU[0][0] - eig1) < epsilon|| fabs(AU[1][1] - eig2) < epsilon)
    {
        cout << "Eigenvalues test passed"<< endl;
    }
    else
    {
        cout << "Eigenvalues test failed"<< endl;
    }
    return;
}

void unity_max()
{
    int n = 5;
    double max_value = 3;
    int k,l;
    double **AU = new double*[n];

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
    if (maxoffdiag == max_value)
    {
        cout << "Maximum non-diagonal test succeded" << endl;
    }
    else
    {
        cout << "Maximum non-diagonal test failed" << endl;
    }
    return;
}

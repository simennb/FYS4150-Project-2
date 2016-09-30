#include "unity.h"
#include "functions.h"
#include <iostream>

void unity(){
    double **AU = new double*[2];
    double **RU = new double*[2];
    for (int i=0; i<n; i++){
        AU[i] = new double[2];
        RU[i] = new double[2];
    }
    AU[1][1] = 2;
    AU[2][2] = 2;
    AU[1][2] = 1;
    AU[2][1] = 1;


}


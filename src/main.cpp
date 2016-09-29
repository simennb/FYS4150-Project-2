#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
    double epsilon = 1.0e-8;

    int iterations = 0;

    while(fabs(maxoffdiag) > epsilon && iterations < maxiterations){
        maxoffdiag = max_offdiag(A, l, k, n);
        rotate();
        iterations++;
    }
}

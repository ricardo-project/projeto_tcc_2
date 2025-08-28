#include<iostream>
#include<iomanip>
#include<cmath>
#include<ctime>
using namespace std;

#include"../radix.hpp"

int main() {

    system("color 4f");
    setlocale(LC_ALL, "Portuguese");
    srand(time(0));

    int i, j, k;
    int Nmx = 1024;
    int N  = 512;
    int N2 = N/2;
    int *I = new int[Nmx];
    Ct *X = new Ct[Nmx];
    Cx *E = new Cx[Nmx];
    setIandE(I, E, Nmx, Nmx/2);

    double v1, v2;
    double *x = new double[N];
    setRand_(x, x + N);

    hartley(x, X, I, E, N, Nmx/N, 0);

    k = 0;
    for(i = 0; i <= N2; i++) {
        v1 = 0;
        v2 = 0;
        for(j = 0; j < N; j++) {
            v1 += x[j]*E[(j*i) % N].a.c;
            v2 += x[j]*E[(j*i) % N].a.s;
        }
        //k += int(abs(v1 - X[i].c) < 1e-10 && abs(v2 - X[i].s) < 1e-10);
    }
    cout << "Quantidade: " << k << " (" << (N2 + 1) << ")\n";
    return 0;
}

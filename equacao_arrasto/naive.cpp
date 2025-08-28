#include<iostream>
#include<iomanip>
using namespace std;

typedef long double ldouble;

int main() {

    int i, j, k;
    int N = 1024;
    ldouble *C  = new ldouble[N + 2];
    ldouble *Cv = new ldouble[N + 1];
    ldouble v1, v2, K = 2;
    cout << setprecision(12);

    C[0] = 1;  // Posição
    C[1] = 1.1;  // Velocidade
    Cv[0] = C[1];

    cout << "C[0] = " << C[0] << endl;
    cout << "C[1] = " << C[1] << endl;

    for(k = 0; k < N; k++) {

        v1 = 0;
        for(j = 0; j <= k; j++) v1 += Cv[j]*Cv[k - j];
        v1 *= K/((k + 1)*(k + 2));

        C[k + 2]  = v1;
        Cv[k + 1] = v1*(k + 1);

        cout << "C[" << (k + 1) << "] = " << v1;
        getchar();
    }
    return 0;
}

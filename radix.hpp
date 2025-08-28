#ifndef HARTLEY
#define HARTLEY

typedef short int sh;
typedef long double ldouble;

struct Ct {
    double c;
    double s;
    double cas;  // cas(x)  = cos(x) + sin(x)
    double casS; // cas*(x) = cos(x) - sin(x)
};

double Rand() {
    return (2*double(rand())/RAND_MAX - 1);
}

double setRand_(double *x, double *xfinal) {
    for(x; x != xfinal; x++) *x = Rand()*0.02;
}

struct Cx {
    Ct a;
};

void setIandE(int *I, Cx *E, int N, int N2) {

    int i, j, k;

    i = N2;
    j = 1;
    I[0] = 0;
    while(j <= N2) {
        for(k = 0; k < j; k++) I[j + k] = I[k] + i;
        i /= 2;
        j *= 2;
    }

    double ang, pi = acos(-1);
    double c, s;
    for(i = 0; i < N; i++) {
        ang = 2*pi*double(i)/N;
        c = cos(ang);
        s = sin(ang);

        (E->a).c = c;
        (E->a).s = s;
        (E->a).cas  = c + s;
        (E->a).casS = c - s;
        E++;
    }
}

void hartley(double *x, Ct *X, int *Ia, Cx *E, int N, int b, bool cond) {

    int i, j, k, l;
    int N4, N2, S, a;
    int *I, *If;
    double Ci, Si, A, B, C, D, F, G, c, s;
    const double r2 = sqrt(2)/2;
    Ct *Xn, *Xa, *XN2, *XN4, *XN2n, *Xfinal = X + N;
    Cx *Eb;
    Cx *Ef = E + (N*b)/4;

    If = Ia + N;
    //for(l = 0; l < Q; l++) {

        Xn = X;
        if(cond) {
            for(I = Ia; I != If; I += 8*b) {

                A = x[ I[4*b] ];
                B = x[ I[6*b] ];
                C = x[ I[2*b] ];
                D = x[ I[0] ];


                // Parte cosseno
                Xn[2].c = D - C;
                s = C + D;
                G = A - B;
                F = G*r2;
                Xn[1].c = D + F;
                Xn[3].c = D - F;

                F = A + B;
                Xn[0].c = s + F;
                Xn[4].c = s - F;


                // Parte seno
                Xn[2].s = G;
                F *= r2;
                Xn[1].s = F + C;
                Xn[3].s = F - C;
                Xn[0].s = 0;
                Xn[4].s = 0;

                Xn += 8;
            }

        } else {
            for(I = Ia; I != If; I += 8) {
                A = x[ I[4*b] ] + x[ I[7*b] ];
                B = x[ I[6*b] ] + x[ I[5*b] ];
                C = x[ I[2*b] ] + x[ I[3*b] ];
                D = x[ I[0] ] + x[ I[b] ];
                F = x[ I[0] ] - x[ I[b] ];

                Xn[2].c = D - C;
                C += D;
                D = (A - B)*r2;
                Xn[1].c = F + D;
                Xn[3].c = F - D;
                A += B;
                Xn[0].c = C + A;
                Xn[4].c = C - A;

                // Parte seno
                A = x[ I[4*b] ] - x[ I[7*b] ];
                B = x[ I[6*b] ] - x[ I[5*b] ];
                C = x[ I[2*b] ] - x[ I[3*b] ];

                Xn[2].s = A - B;
                A = (A + B)*r2;
                Xn[1].s = A + C;
                Xn[3].s = A - C;
                Xn[0].s = 0;
                Xn[4].s = 0;
                Xn += 8;
            }
        }

        N4 = 4;
        N2 = 2*N4;
        S = 2*N2;
        a = (N*b)/S;

        while(S <= N) {
            for(Xn = X; Xn != Xfinal; Xn += S) {
                Xa = Xn;
                XN2  = Xn + N2;
                XN2n = XN2;

                A = XN2->c;
                B = XN2->s;

                XN2n->c =  Xa->c - A;
                XN2n->s = -Xa->s + B;
                Xa->c += A;
                Xa->s += B;

                Xa++;
                XN2++;
                XN2n--;

                for(Eb = E + a; Eb != Ef; Eb += a) {
                    Ci = Xa->c;
                    Si = Xa->s;

                    c = Eb->a.c;
                    s = Eb->a.s;

                    A = XN2->c*c - XN2->s*s;
                    B = XN2->c*s + XN2->s*c;

                    Xa->c = Ci + A;
                    Xa->s = Si + B;

                    XN2n->c =  Ci - A;
                    XN2n->s = -Si + B;

                    Xa++;
                    XN2++;
                    XN2n--;
                }

                // FINAL DO PROCESSO
                // *****************

                Ci = Xa->c;
                Si = Xa->s;

                A = -XN2->s;
                B =  XN2->c;

                XN2n->c =  Ci - A;
                XN2n->s = -Si + B;
                Xa->c = Ci + A;
                Xa->s = Si + B;
            }
            N4 *= 2;
            N2 *= 2;
            S *= 2;
            a /= 2;
        }
    //}
}

void Hadamard(double *z, Ct *W, Ct *Z, double *zfinal, int N) {
    *z = W->c*Z->c;
    W++;
    Z++;
    double *zi;
    double *zNi = z + N - 1;
    double A, B;
    for(zi = z + 1; zi != zfinal; zi++) {
        A = W->c*Z->c - W->s*Z->s;
        B = W->c*Z->s + W->s*Z->c;
        *zi  = A + B;
        *zNi = A - B;
        zNi--;
        Z++;
        W++;
    } *zi = W->c*Z->c;
}

void HadamardCR(double *z, Ct *W, Ct *Z, double *zfinal, int N) { // Cross-relation
    *z = W->c*Z->c;
    W++;
    Z++;
    double *zi;
    double *zNi = z + N - 1;
    double A, B;
    for(zi = z + 1; zi != zfinal; zi++) {
        A = W->c*Z->c + W->s*Z->s;
        B = W->c*Z->s - W->s*Z->c;
        *zi  = A + B;
        *zNi = A - B;
        zNi--;
        Z++;
        W++;
    } *zi = W->c*Z->c;
}

void convFinal(double *z, Ct *Z, double *zfinal, int N) {
    *z = Z->c/N;
    Ct *Zi = Z + 1;
    double *zi;
    double *zNi = z + N - 1;
    for(zi = z + 1; zi != zfinal; zi++) {
        *zi = (Zi->c + Zi->s)/N;
        *zNi = (Zi->c - Zi->s)/N;
        Zi++;
        zNi--;
    } *zi = Zi->c/N;
}

void addHartley(double *a, Ct *W, int d2) {
    a[0]  += W[0].c;
    a[d2] += W[d2].c;
    for(int k = 1; k < d2; k++) a[k] += W[k].c + W[k].s;
}


#endif // HARTLEY

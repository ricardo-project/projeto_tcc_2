#ifndef HARTLEY
#define HARTLEY

typedef short int sh;

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

void hartley(double *x, Ct *X, int *Ia, Cx *E, int N, int Q, bool cond) {

    int i, j, k, l;
    int N4, N2, S, a;
    int *I, *If;
    double Ci, Si, A, B, C, D, F, G, c, s;
    const double r2 = sqrt(2)/2;
    Ct *Xn, *Xa, *XN2, *XN4, *XN2n, *Xfinal = X + N;
    Cx *Eb;
    Cx *Ef = E + N/4;
    Cx *Eg = E + N/8;

    If = Ia + N;
    for(l = 0; l < Q; l++) {

        Xn = X;
        if(cond) {
            for(I = Ia; I != If; I += 8) {

                A = x[ I[4] ];
                B = x[ I[6] ];
                C = x[ I[2] ];
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
                A = x[ I[4] ] + x[ I[7] ];
                B = x[ I[6] ] + x[ I[5] ];
                C = x[ I[2] ] + x[ I[3] ];
                D = x[ I[0] ] + x[ I[1] ];
                F = x[ I[0] ] - x[ I[1] ];

                Xn[2].c = D - C;
                C += D;
                D = (A - B)*r2;
                Xn[1].c = F + D;
                Xn[3].c = F - D;
                A += B;
                Xn[0].c = C + A;
                Xn[4].c = C - A;

                // Parte seno
                A = x[ I[4] ] - x[ I[7] ];
                B = x[ I[6] ] - x[ I[5] ];
                C = x[ I[2] ] - x[ I[3] ];

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
        a = N/S;

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
    }
}

void halfHartley(double *x, Ct *X, int *Ia, Cx *E, int N, int Q) {

    int i, j, k, l;
    int N4, N2, S, a;
    int N2c = N/2;
    int *I, *If;
    double Ci, Si, A, B, C, D, F, c, s;
    double r2 = sqrt(2)/2;
    Ct *Xn, *Xa, *XN2, *XN4, *XN2n, *Xfinal = X + N2c;
    Cx *Eb;
    Cx *Ef = E + N/4;
    Cx *Eg;

    If = Ia + N;
    for(l = 0; l < Q; l++) {

        Xn = X;
        for(I = Ia; I != If; I += 8) {
            A = x[  *I  ] - x[ I[1] ];
            B = x[ I[2] ] - x[ I[3] ];
            C = x[ I[4] ] - x[ I[5] ];
            F = x[ I[6] ] - x[ I[7] ];

            // Parte cosseno
            D = (C - F)*r2;
            Xn[0].c = A + D;
            Xn[1].c = A - D;

            // Parte seno
            A = (C + F)*r2;
            Xn[0].s = A + B;
            Xn[1].s = A - B;
            Xn += 4;
        }

        N4 = 2;
        N2 = 2*N4;
        S = 2*N2;
        a = N/S;

        while(S <= N2c) {
            Eg = Ef + a/2;
            for(Xn = X; Xn != Xfinal; Xn += S) {
                Xa = Xn;
                XN2  = Xn + N2;
                XN2n = XN2 - 1;

                for(Eb = E + a/2; Eb != Eg; Eb += a) {
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
            }
            N4 *= 2;
            N2 *= 2;
            S *= 2;
            a /= 2;
        }
    }
}

void halfHadamard(double *z, Ct *W, Ct *Z, double *zfinal, int N2) {
    int k;
    double V, L;
    double *zNi = z + N2 - 1;
    for(z; z != zfinal; z++) {
        V = W->c + W->s;
        L = W->c - W->s;
        *z   = Z->c*V + Z->s*L;
        *zNi = Z->c*L - Z->s*V;
        zNi--;
        W++;
        Z++;
    }
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

void halfConvFinal(double *z, double *g, Ct *X, Ct *Xfinal, Cx *E, int N, int N2, int N4) {
    int k;
    double V, L;

    Ct *Xini = X;
    *z = *g - 2*X->c/N;
    g++;
    z++;
    E++;
    for(X++; X != Xfinal; X++) {
        *z = *g - 2*(E->a.cas*X->c + E->a.casS*X->s)/N;
        g++;
        z++;
        E++;
    }
    *z = *g - 2*X->c/N;
    g++;
    z++;
    E++;

    for(X--; X != Xini; X--) {
        *z = *g - 2*(E->a.cas*X->c - E->a.casS*X->s)/N;
        g++;
        z++;
        E++;
    }
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

void halfConvEnd(double *z, double *g, Ct *Z, double *zfinal, int N) {
    //*z = Z->c/N;

    *g = Z->c/N - (*z);
    Ct *Zi = Z + 1;
    double *zi;
    double *gi = g + 1;
    double *zNi = z + N - 1;
    double *gNi = g + N - 1;

    for(zi = z + 1; zi != zfinal; zi++) {
        *gi = (Zi->c + Zi->s)/N - (*zi);
        *gNi = (Zi->c - Zi->s)/N - (*zNi);
        Zi++;
        gi++;
        zNi--;
        gNi--;
    } *gi = Zi->c/N - (*zi);
}


#endif // HARTLEY

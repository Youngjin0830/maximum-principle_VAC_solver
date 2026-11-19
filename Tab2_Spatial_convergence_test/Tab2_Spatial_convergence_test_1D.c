#pragma STDC FP_CONTRACT OFF
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ID(i,p) ((size_t)(i)+Sx*(size_t)(p))

int main(void)
{
    const double pi = 3.141592653589793;
    const double Lx = 0, Rx = 1;
    const int N = 4;
    const double eps = 0.05;
    const double T = 1.0e-3;
    const int MM[6] = {32, 64, 128, 256, 512, 1024};
    const int nM = 6;

    double *sol[6];
    int nD = 0;
    for (int l = 0; l < nM; l++) {
        const int Nx = MM[l]; const double h = (Rx-Lx)/Nx;
        const size_t Sx = Nx+2;
        const size_t V = Sx*N;
        double x[Nx+2];
        for (int i = 0; i < Nx+2; i++) x[i] = (Lx-0.5*h)+(i*((Rx+0.5*h)-(Lx-0.5*h)))/(Nx+1);
        x[0] = Lx-0.5*h; x[Nx+1] = Rx+0.5*h;
        const int Nt = 20*(Nx/32)*(Nx/32); const double dt = T/Nt;
        const double e2 = eps*eps, h2 = h*h;
        printf("1D M = %d: dt = %.15e, dtmax = %.15e\n", Nx, dt, 2*e2*h2/(h2+4*e2));
        fflush(stdout);

        double *u = malloc(V*sizeof(double)), *w = malloc(V*sizeof(double)), *t;
        sol[l] = malloc((size_t)Nx*N*sizeof(double));
        if (!u || !w || !sol[l]) { printf("1D M = %d: not enough memory\n", Nx); free(u); free(w); free(sol[l]); break; }

        for (int i = 0; i < Nx+2; i++) {
            u[ID(i,0)] = 0.25+0.1*cos(pi*x[i]); u[ID(i,1)] = 0.25+0.1*cos(2*pi*x[i]);
            u[ID(i,2)] = 0.25-0.1*cos(pi*x[i]); u[ID(i,3)] = 0.25-0.1*cos(2*pi*x[i]);
        }
        for (int p = 0; p < N; p++) { u[ID(0,p)] = u[ID(1,p)]; u[ID(Nx+1,p)] = u[ID(Nx,p)]; }
        for (size_t q = 0; q < V; q++) w[q] = u[q];
        for (int it = 1; it <= Nt; it++) {
            for (int i = 1; i <= Nx; i++) {
                double Fp[4], beta;
                for (int p = 0; p < N; p++) { double v = u[ID(i,p)]; Fp[p] = v*(v-0.5)*(v-1); }
                beta = Fp[0];
                for (int p = 1; p < N; p++) beta = beta+Fp[p];
                beta = beta/(N*e2);
                for (int p = 0; p < N; p++)
                    w[ID(i,p)] = u[ID(i,p)]+dt*(-Fp[p]/e2+(u[ID(i-1,p)]+u[ID(i+1,p)]-2*u[ID(i,p)])/h2+beta);
            }
            for (int p = 0; p < N; p++) { w[ID(0,p)] = w[ID(1,p)]; w[ID(Nx+1,p)] = w[ID(Nx,p)]; }
            t = u; u = w; w = t;
        }
        for (int p = 0; p < N; p++) for (int i = 1; i <= Nx; i++) sol[l][(i-1)+(size_t)Nx*p] = u[ID(i,p)];
        free(u); free(w);
        nD = l+1;
    }

    double err[5], rate[5];
    for (int l = 0; l < nD-1; l++) {
        const double h = (Rx-Lx)/MM[l];
        const int Mc = MM[l], Mf = MM[l+1];
        const double *pf = sol[l+1];
        double sum = 0;
        for (int p = 0; p < N; p++) for (int a = 0; a < Mc; a++) {
            double pc = (pf[2*a+(size_t)Mf*p]+pf[2*a+1+(size_t)Mf*p])/2;
            double d = sol[l][a+(size_t)Mc*p]-pc;
            sum = sum+d*d;
        }
        err[l] = sqrt(h*sum/N);
        printf("1D M = %d: err = %.15e\n", MM[l], err[l]);
    }
    for (int l = 1; l < nD-1; l++) rate[l] = log2(err[l-1]/err[l]);
    printf("1D M    = ");
    for (int l = 0; l < nD-1; l++) printf("%d ", MM[l]);
    printf("\n1D err  = ");
    for (int l = 0; l < nD-1; l++) printf("%.4e ", err[l]);
    printf("\n1D rate = ");
    for (int l = 0; l < nD-1; l++) { if (l == 0) printf("NaN "); else printf("%.4f ", rate[l]); }
    printf("\n");

    for (int l = 0; l < nD; l++) free(sol[l]);
    return 0;
}

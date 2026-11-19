#pragma STDC FP_CONTRACT OFF
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ID(i,j,p) ((size_t)(i)+Sx*((size_t)(j)+Sy*(size_t)(p)))

int main(void)
{
    const double pi = 3.141592653589793;
    const double Lx = 0, Rx = 1, Ly = 0, Ry = 1;
    const int N = 4;
    const double eps = 0.05;
    const double T = 1.0e-3;
    const int MM[6] = {32, 64, 128, 256, 512, 1024};
    const int nM = 6;

    double *sol[6];
    int nD = 0;
    for (int l = 0; l < nM; l++) {
        const int Nx = MM[l]; const double h = (Rx-Lx)/Nx; const int Ny = (int)((Ry-Ly)/h);
        const size_t Sx = Nx+2, Sy = Ny+2;
        const size_t V = Sx*Sy*N;
        double x[Nx+2], y[Ny+2];
        for (int i = 0; i < Nx+2; i++) x[i] = (Lx-0.5*h)+(i*((Rx+0.5*h)-(Lx-0.5*h)))/(Nx+1);
        for (int j = 0; j < Ny+2; j++) y[j] = (Ly-0.5*h)+(j*((Ry+0.5*h)-(Ly-0.5*h)))/(Ny+1);
        x[0] = Lx-0.5*h; x[Nx+1] = Rx+0.5*h;
        y[0] = Ly-0.5*h; y[Ny+1] = Ry+0.5*h;
        const int Nt = 20*(Nx/32)*(Nx/32); const double dt = T/Nt;
        const double e2 = eps*eps, h2 = h*h;
        printf("2D M = %d: dt = %.15e, dtmax = %.15e\n", Nx, dt, 2*e2*h2/(h2+8*e2));
        fflush(stdout);

        double *u = malloc(V*sizeof(double)), *w = malloc(V*sizeof(double)), *t;
        sol[l] = malloc((size_t)Nx*Ny*N*sizeof(double));
        if (!u || !w || !sol[l]) { printf("2D M = %d: not enough memory\n", Nx); free(u); free(w); free(sol[l]); break; }

        for (int j = 0; j < Ny+2; j++) for (int i = 0; i < Nx+2; i++) {
            double c1 = cos(pi*x[i])*cos(pi*y[j]), c2 = cos(2*pi*x[i])*cos(2*pi*y[j]);
            u[ID(i,j,0)] = 0.25+0.1*c1; u[ID(i,j,1)] = 0.25+0.1*c2;
            u[ID(i,j,2)] = 0.25-0.1*c1; u[ID(i,j,3)] = 0.25-0.1*c2;
        }
        for (int p = 0; p < N; p++) for (int j = 0; j < Ny+2; j++) { u[ID(0,j,p)] = u[ID(1,j,p)]; u[ID(Nx+1,j,p)] = u[ID(Nx,j,p)]; }
        for (int p = 0; p < N; p++) for (int i = 0; i < Nx+2; i++) { u[ID(i,0,p)] = u[ID(i,1,p)]; u[ID(i,Ny+1,p)] = u[ID(i,Ny,p)]; }
        for (size_t q = 0; q < V; q++) w[q] = u[q];
        for (int it = 1; it <= Nt; it++) {
            #pragma omp parallel for
            for (int j = 1; j <= Ny; j++) for (int i = 1; i <= Nx; i++) {
                double Fp[4], beta;
                for (int p = 0; p < N; p++) { double v = u[ID(i,j,p)]; Fp[p] = v*(v-0.5)*(v-1); }
                beta = Fp[0];
                for (int p = 1; p < N; p++) beta = beta+Fp[p];
                beta = beta/(N*e2);
                for (int p = 0; p < N; p++)
                    w[ID(i,j,p)] = u[ID(i,j,p)]+dt*(-Fp[p]/e2+(u[ID(i-1,j,p)]+u[ID(i+1,j,p)]+u[ID(i,j-1,p)]+u[ID(i,j+1,p)]-4*u[ID(i,j,p)])/h2+beta);
            }
            for (int p = 0; p < N; p++) for (int j = 0; j < Ny+2; j++) { w[ID(0,j,p)] = w[ID(1,j,p)]; w[ID(Nx+1,j,p)] = w[ID(Nx,j,p)]; }
            for (int p = 0; p < N; p++) for (int i = 0; i < Nx+2; i++) { w[ID(i,0,p)] = w[ID(i,1,p)]; w[ID(i,Ny+1,p)] = w[ID(i,Ny,p)]; }
            t = u; u = w; w = t;
        }
        for (int p = 0; p < N; p++) for (int j = 1; j <= Ny; j++) for (int i = 1; i <= Nx; i++)
            sol[l][(i-1)+(size_t)Nx*((j-1)+(size_t)Ny*p)] = u[ID(i,j,p)];
        free(u); free(w);
        nD = l+1;
    }

    double err[5], rate[5];
    for (int l = 0; l < nD-1; l++) {
        const double h = (Rx-Lx)/MM[l];
        const int Mc = MM[l], Mf = MM[l+1];
        const double *pf = sol[l+1];
        double sum = 0;
        for (int p = 0; p < N; p++) for (int b = 0; b < Mc; b++) for (int a = 0; a < Mc; a++) {
            const size_t f = (size_t)Mf*Mf*p;
            double pc = (pf[f+2*a+(size_t)Mf*(2*b)]+pf[f+2*a+1+(size_t)Mf*(2*b)]+pf[f+2*a+(size_t)Mf*(2*b+1)]+pf[f+2*a+1+(size_t)Mf*(2*b+1)])/4;
            double d = sol[l][a+(size_t)Mc*(b+(size_t)Mc*p)]-pc;
            sum = sum+d*d;
        }
        err[l] = sqrt(h*h*sum/N);
        printf("2D M = %d: err = %.15e\n", MM[l], err[l]);
    }
    for (int l = 1; l < nD-1; l++) rate[l] = log2(err[l-1]/err[l]);
    printf("2D M    = ");
    for (int l = 0; l < nD-1; l++) printf("%d ", MM[l]);
    printf("\n2D err  = ");
    for (int l = 0; l < nD-1; l++) printf("%.4e ", err[l]);
    printf("\n2D rate = ");
    for (int l = 0; l < nD-1; l++) { if (l == 0) printf("NaN "); else printf("%.4f ", rate[l]); }
    printf("\n");

    for (int l = 0; l < nD; l++) free(sol[l]);
    return 0;
}

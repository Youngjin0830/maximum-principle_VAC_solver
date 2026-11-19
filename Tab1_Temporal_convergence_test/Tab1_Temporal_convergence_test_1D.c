#pragma STDC FP_CONTRACT OFF
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ID(i,p) ((size_t)(i)+Sx*(size_t)(p))

int main(void)
{
    const double pi = 3.141592653589793;
    const int Nx = 256;
    const double Lx = 0, Rx = 1, h = (Rx-Lx)/Nx;
    const int N = 4;
    const size_t Sx = Nx+2;
    const size_t V = Sx*N;
    double x[Nx+2];
    for (int i = 0; i < Nx+2; i++) x[i] = (Lx-0.5*h)+(i*((Rx+0.5*h)-(Lx-0.5*h)))/(Nx+1);
    x[0] = Lx-0.5*h; x[Nx+1] = Rx+0.5*h;

    const double m = 4;
    const double eps = m*h/(4*sqrt(2.0)*atanh(0.9));
    const double e2 = eps*eps, h2 = h*h;

    const double dtmax = 2*e2*h2/(h2+4*e2);

    double *phi0 = malloc(V*sizeof(double));
    double *A = malloc(V*sizeof(double)), *B = malloc(V*sizeof(double));
    double *P = malloc(V*sizeof(double)), *Q = malloc(V*sizeof(double));
    if (!phi0 || !A || !B || !P || !Q) { printf("1D: not enough memory\n"); return 1; }

    for (int i = 0; i < Nx+2; i++) {
        double c1 = cos(pi*x[i]), c2 = cos(2*pi*x[i]);
        phi0[ID(i,0)] = 0.25+0.1*c1; phi0[ID(i,1)] = 0.25+0.1*c2;
        phi0[ID(i,2)] = 0.25-0.1*c1; phi0[ID(i,3)] = 0.25-0.1*c2;
    }
    for (int p = 0; p < N; p++) { phi0[ID(0,p)] = phi0[ID(1,p)]; phi0[ID(Nx+1,p)] = phi0[ID(Nx,p)]; }

    double dt = dtmax/pow(2,14);
    for (size_t q = 0; q < V; q++) { A[q] = phi0[q]; B[q] = phi0[q]; }
    double *u = A, *w = B, *t;
    for (int it = 1; it <= 16*(1<<14); it++) {
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
    double *ref = u;

    const int nL = 6;
    double err[6], dts[6], rate[5];
    for (int il = 0; il < nL; il++) {
        const int L = il+1;
        dt = dtmax/pow(2,L); dts[il] = dt;
        for (size_t q = 0; q < V; q++) { P[q] = phi0[q]; Q[q] = phi0[q]; }
        u = P; w = Q;
        for (int it = 1; it <= 16*(1<<L); it++) {
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
        double sum = 0;
        for (int p = 0; p < N; p++) for (int i = 1; i <= Nx; i++) {
            double d = ref[ID(i,p)]-u[ID(i,p)];
            sum = sum+d*d;
        }
        err[il] = sqrt(sum/((double)N*Nx));
        printf("1D dt = %.15e: err = %.15e\n", dts[il], err[il]);
        fflush(stdout);
    }
    for (int il = 0; il < nL-1; il++) rate[il] = log2(err[il]/err[il+1]);
    printf("1D err  = ");
    for (int il = 0; il < nL; il++) printf("%.4e ", err[il]);
    printf("\n1D rate = ");
    for (int il = 0; il < nL-1; il++) printf("%.4f ", rate[il]);
    printf("\n");

    free(phi0); free(A); free(B); free(P); free(Q);
    return 0;
}

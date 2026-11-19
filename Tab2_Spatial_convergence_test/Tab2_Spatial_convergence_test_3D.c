#pragma STDC FP_CONTRACT OFF
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ID(i,j,k,p) ((size_t)(i)+Sx*((size_t)(j)+Sy*((size_t)(k)+Sz*(size_t)(p))))

int main(void)
{
    const double pi = 3.141592653589793;
    const double Lx = 0, Rx = 1, Ly = 0, Ry = 1, Lz = 0, Rz = 1;
    const int N = 4;
    const double eps = 0.05;
    const double T = 1.0e-3;
    const int MM[6] = {32, 64, 128, 256, 512, 1024};
    const int nM = 6;

    double err[5], rate[5];
    double *prev = NULL;
    int nE = 0;
    for (int l = 0; l < nM; l++) {
        const int Nx = MM[l]; const double h = (Rx-Lx)/Nx; const int Ny = (int)((Ry-Ly)/h), Nz = (int)((Rz-Lz)/h);
        const size_t Sx = Nx+2, Sy = Ny+2, Sz = Nz+2;
        const size_t V = Sx*Sy*Sz*N;
        double x[Nx+2], y[Ny+2], z[Nz+2];
        for (int i = 0; i < Nx+2; i++) x[i] = (Lx-0.5*h)+(i*((Rx+0.5*h)-(Lx-0.5*h)))/(Nx+1);
        for (int j = 0; j < Ny+2; j++) y[j] = (Ly-0.5*h)+(j*((Ry+0.5*h)-(Ly-0.5*h)))/(Ny+1);
        for (int k = 0; k < Nz+2; k++) z[k] = (Lz-0.5*h)+(k*((Rz+0.5*h)-(Lz-0.5*h)))/(Nz+1);
        x[0] = Lx-0.5*h; x[Nx+1] = Rx+0.5*h;
        y[0] = Ly-0.5*h; y[Ny+1] = Ry+0.5*h;
        z[0] = Lz-0.5*h; z[Nz+1] = Rz+0.5*h;
        const int Nt = 20*(Nx/32)*(Nx/32); const double dt = T/Nt;
        const double e2 = eps*eps, h2 = h*h;
        printf("3D M = %d: dt = %.15e, dtmax = %.15e\n", Nx, dt, 2*e2*h2/(h2+12*e2));
        fflush(stdout);

        double *u = malloc(V*sizeof(double)), *w = malloc(V*sizeof(double)), *t;
        if (!u || !w) { printf("3D M = %d: not enough memory (%.1f GB per array)\n", Nx, V*sizeof(double)/1e9); free(u); free(w); break; }

        for (int k = 0; k < Nz+2; k++) for (int j = 0; j < Ny+2; j++) for (int i = 0; i < Nx+2; i++) {
            double c1 = cos(pi*x[i])*cos(pi*y[j])*cos(pi*z[k]), c2 = cos(2*pi*x[i])*cos(2*pi*y[j])*cos(2*pi*z[k]);
            u[ID(i,j,k,0)] = 0.25+0.1*c1; u[ID(i,j,k,1)] = 0.25+0.1*c2;
            u[ID(i,j,k,2)] = 0.25-0.1*c1; u[ID(i,j,k,3)] = 0.25-0.1*c2;
        }
        for (int p = 0; p < N; p++) for (int k = 0; k < Nz+2; k++) for (int j = 0; j < Ny+2; j++) { u[ID(0,j,k,p)] = u[ID(1,j,k,p)]; u[ID(Nx+1,j,k,p)] = u[ID(Nx,j,k,p)]; }
        for (int p = 0; p < N; p++) for (int k = 0; k < Nz+2; k++) for (int i = 0; i < Nx+2; i++) { u[ID(i,0,k,p)] = u[ID(i,1,k,p)]; u[ID(i,Ny+1,k,p)] = u[ID(i,Ny,k,p)]; }
        for (int p = 0; p < N; p++) for (int j = 0; j < Ny+2; j++) for (int i = 0; i < Nx+2; i++) { u[ID(i,j,0,p)] = u[ID(i,j,1,p)]; u[ID(i,j,Nz+1,p)] = u[ID(i,j,Nz,p)]; }
        for (size_t q = 0; q < V; q++) w[q] = u[q];
        for (int it = 1; it <= Nt; it++) {
            #pragma omp parallel for
            for (int k = 1; k <= Nz; k++) for (int j = 1; j <= Ny; j++) for (int i = 1; i <= Nx; i++) {
                double Fp[4], beta;
                for (int p = 0; p < N; p++) { double v = u[ID(i,j,k,p)]; Fp[p] = v*(v-0.5)*(v-1); }
                beta = Fp[0];
                for (int p = 1; p < N; p++) beta = beta+Fp[p];
                beta = beta/(N*e2);
                for (int p = 0; p < N; p++)
                    w[ID(i,j,k,p)] = u[ID(i,j,k,p)]+dt*(-Fp[p]/e2+(u[ID(i-1,j,k,p)]+u[ID(i+1,j,k,p)]+u[ID(i,j-1,k,p)]+u[ID(i,j+1,k,p)]+u[ID(i,j,k-1,p)]+u[ID(i,j,k+1,p)]-6*u[ID(i,j,k,p)])/h2+beta);
            }
            for (int p = 0; p < N; p++) for (int k = 0; k < Nz+2; k++) for (int j = 0; j < Ny+2; j++) { w[ID(0,j,k,p)] = w[ID(1,j,k,p)]; w[ID(Nx+1,j,k,p)] = w[ID(Nx,j,k,p)]; }
            for (int p = 0; p < N; p++) for (int k = 0; k < Nz+2; k++) for (int i = 0; i < Nx+2; i++) { w[ID(i,0,k,p)] = w[ID(i,1,k,p)]; w[ID(i,Ny+1,k,p)] = w[ID(i,Ny,k,p)]; }
            for (int p = 0; p < N; p++) for (int j = 0; j < Ny+2; j++) for (int i = 0; i < Nx+2; i++) { w[ID(i,j,0,p)] = w[ID(i,j,1,p)]; w[ID(i,j,Nz+1,p)] = w[ID(i,j,Nz,p)]; }
            t = u; u = w; w = t;
        }
        free(w);

        if (l > 0) {
            const int Mc = MM[l-1];
            double sum = 0;
            for (int p = 0; p < N; p++) for (int c = 0; c < Mc; c++) for (int b = 0; b < Mc; b++) for (int a = 0; a < Mc; a++) {
                const int i1 = 2*a+1, i2 = 2*a+2, j1 = 2*b+1, j2 = 2*b+2, k1 = 2*c+1, k2 = 2*c+2;
                double pc = (u[ID(i1,j1,k1,p)]+u[ID(i2,j1,k1,p)]+u[ID(i1,j2,k1,p)]+u[ID(i2,j2,k1,p)]
                            +u[ID(i1,j1,k2,p)]+u[ID(i2,j1,k2,p)]+u[ID(i1,j2,k2,p)]+u[ID(i2,j2,k2,p)])/8;
                double d = prev[a+(size_t)Mc*(b+(size_t)Mc*(c+(size_t)Mc*p))]-pc;
                sum = sum+d*d;
            }
            err[l-1] = sqrt(pow(2*h,3)*sum/N);
            nE = l;
            printf("3D M = %d: err = %.15e\n", Mc, err[l-1]);
            fflush(stdout);
            free(prev); prev = NULL;
        }
        if (l < nM-1) {
            prev = malloc((size_t)Nx*Ny*Nz*N*sizeof(double));
            if (!prev) { printf("3D M = %d: not enough memory\n", Nx); free(u); break; }
            for (int p = 0; p < N; p++) for (int k = 1; k <= Nz; k++) for (int j = 1; j <= Ny; j++) for (int i = 1; i <= Nx; i++)
                prev[(i-1)+(size_t)Nx*((j-1)+(size_t)Ny*((k-1)+(size_t)Nz*p))] = u[ID(i,j,k,p)];
        }
        free(u);
    }
    free(prev);

    for (int l = 1; l < nE; l++) rate[l] = log2(err[l-1]/err[l]);
    printf("3D M    = ");
    for (int l = 0; l < nE; l++) printf("%d ", MM[l]);
    printf("\n3D err  = ");
    for (int l = 0; l < nE; l++) printf("%.4e ", err[l]);
    printf("\n3D rate = ");
    for (int l = 0; l < nE; l++) { if (l == 0) printf("NaN "); else printf("%.4f ", rate[l]); }
    printf("\n");
    return 0;
}

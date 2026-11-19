clc; clear all; close all;

Nx = 256;
Lx = 0; Rx = 1; h = (Rx-Lx)/Nx;
x = linspace(Lx-0.5*h,Rx+0.5*h,Nx+2)';
N = 4;
I = 2:Nx+1;

m = 4;
eps = m*h/(4*sqrt(2)*atanh(0.9));

dtmax = 2*eps^2*h^2/(h^2+4*eps^2);

c1 = cos(pi*x); c2 = cos(2*pi*x);
phi0 = [1/4+0.1*c1,1/4+0.1*c2,1/4-0.1*c1,1/4-0.1*c2];
phi0(1,:) = phi0(2,:); phi0(Nx+2,:) = phi0(Nx+1,:);

dt = dtmax/2^14;
ref = phi0;
for it = 1:16*2^14
    Fp = ref(I,:).*(ref(I,:)-0.5).*(ref(I,:)-1);
    beta = sum(Fp,2)/(N*eps^2);
    ref(I,:) = ref(I,:)+dt*(-Fp/eps^2+(ref(1:Nx,:)+ref(3:Nx+2,:)-2*ref(I,:))/h^2+beta);
    ref(1,:) = ref(2,:); ref(Nx+2,:) = ref(Nx+1,:);
end

L = 1:6;
err = zeros(1,numel(L)); dts = err;
for il = 1:numel(L)
dt = dtmax/2^L(il); dts(il) = dt;
phi = phi0;
for it = 1:16*2^L(il)
    Fp = phi(I,:).*(phi(I,:)-0.5).*(phi(I,:)-1);
    beta = sum(Fp,2)/(N*eps^2);
    phi(I,:) = phi(I,:)+dt*(-Fp/eps^2+(phi(1:Nx,:)+phi(3:Nx+2,:)-2*phi(I,:))/h^2+beta);
    phi(1,:) = phi(2,:); phi(Nx+2,:) = phi(Nx+1,:);
end
err(il) = sqrt(sum((ref(I,:)-phi(I,:)).^2,'all')/(N*Nx));
end
rate = log2(err(1:end-1)./err(2:end));
fprintf('1D err  = %s\n',sprintf('%.4e ',err));
fprintf('1D rate = %s\n',sprintf('%.4f ',rate));

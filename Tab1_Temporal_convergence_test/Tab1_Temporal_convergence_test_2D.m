clc; clear all; close all;

Nx = 256;
Lx = 0; Rx = 1; h = (Rx-Lx)/Nx;
Ly = 0; Ry = 1; Ny = (Ry-Ly)/h;
x = linspace(Lx-0.5*h,Rx+0.5*h,Nx+2);
y = linspace(Ly-0.5*h,Ry+0.5*h,Ny+2);
N = 4;
I = 2:Nx+1; J = 2:Ny+1;

m = 4;
eps = m*h/(4*sqrt(2)*atanh(0.9));

dtmax = 2*eps^2*h^2/(h^2+8*eps^2);

[xx,yy] = ndgrid(x,y);
c1 = cos(pi*xx).*cos(pi*yy); c2 = cos(2*pi*xx).*cos(2*pi*yy);
phi0 = cat(3,1/4+0.1*c1,1/4+0.1*c2,1/4-0.1*c1,1/4-0.1*c2);
phi0(1,:,:) = phi0(2,:,:); phi0(Nx+2,:,:) = phi0(Nx+1,:,:);
phi0(:,1,:) = phi0(:,2,:); phi0(:,Ny+2,:) = phi0(:,Ny+1,:);

dt = dtmax/2^14;
ref = phi0;
for it = 1:16*2^14
    Fp = ref(I,J,:).*(ref(I,J,:)-0.5).*(ref(I,J,:)-1);
    beta = sum(Fp,3)/(N*eps^2);
    ref(I,J,:) = ref(I,J,:)+dt*(-Fp/eps^2+(ref(1:Nx,J,:)+ref(3:Nx+2,J,:)+ref(I,1:Ny,:)+ref(I,3:Ny+2,:)-4*ref(I,J,:))/h^2+beta);
    ref(1,:,:) = ref(2,:,:); ref(Nx+2,:,:) = ref(Nx+1,:,:);
    ref(:,1,:) = ref(:,2,:); ref(:,Ny+2,:) = ref(:,Ny+1,:);
end

L = 1:6;
err = zeros(1,numel(L)); dts = err;
for il = 1:numel(L)
dt = dtmax/2^L(il); dts(il) = dt;
phi = phi0;
for it = 1:16*2^L(il)
    Fp = phi(I,J,:).*(phi(I,J,:)-0.5).*(phi(I,J,:)-1);
    beta = sum(Fp,3)/(N*eps^2);
    phi(I,J,:) = phi(I,J,:)+dt*(-Fp/eps^2+(phi(1:Nx,J,:)+phi(3:Nx+2,J,:)+phi(I,1:Ny,:)+phi(I,3:Ny+2,:)-4*phi(I,J,:))/h^2+beta);
    phi(1,:,:) = phi(2,:,:); phi(Nx+2,:,:) = phi(Nx+1,:,:);
    phi(:,1,:) = phi(:,2,:); phi(:,Ny+2,:) = phi(:,Ny+1,:);
end
err(il) = sqrt(sum((ref(I,J,:)-phi(I,J,:)).^2,'all')/(N*Nx*Ny));
end
rate = log2(err(1:end-1)./err(2:end));
fprintf('2D err  = %s\n',sprintf('%.4e ',err));
fprintf('2D rate = %s\n',sprintf('%.4f ',rate));

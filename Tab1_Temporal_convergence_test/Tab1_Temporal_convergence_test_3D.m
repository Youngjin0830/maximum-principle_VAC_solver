clc; clear all; close all;

Nx = 256;
Lx = 0; Rx = 1; h = (Rx-Lx)/Nx;
Ly = 0; Ry = 1; Ny = (Ry-Ly)/h;
Lz = 0; Rz = 1; Nz = (Rz-Lz)/h;
x = linspace(Lx-0.5*h,Rx+0.5*h,Nx+2);
y = linspace(Ly-0.5*h,Ry+0.5*h,Ny+2);
z = linspace(Lz-0.5*h,Rz+0.5*h,Nz+2);
N = 4;
I = 2:Nx+1; J = 2:Ny+1; K = 2:Nz+1;

m = 4;
eps = m*h/(4*sqrt(2)*atanh(0.9));

dtmax = 2*eps^2*h^2/(h^2+12*eps^2);

[xx,yy,zz] = ndgrid(x,y,z);
c1 = cos(pi*xx).*cos(pi*yy).*cos(pi*zz); c2 = cos(2*pi*xx).*cos(2*pi*yy).*cos(2*pi*zz);
phi0 = cat(4,1/4+0.1*c1,1/4+0.1*c2,1/4-0.1*c1,1/4-0.1*c2);
phi0(1,:,:,:) = phi0(2,:,:,:); phi0(Nx+2,:,:,:) = phi0(Nx+1,:,:,:);
phi0(:,1,:,:) = phi0(:,2,:,:); phi0(:,Ny+2,:,:) = phi0(:,Ny+1,:,:);
phi0(:,:,1,:) = phi0(:,:,2,:); phi0(:,:,Nz+2,:) = phi0(:,:,Nz+1,:);

dt = dtmax/2^14;
ref = phi0;
for it = 1:16*2^14
    Fp = ref(I,J,K,:).*(ref(I,J,K,:)-0.5).*(ref(I,J,K,:)-1);
    beta = sum(Fp,4)/(N*eps^2);
    ref(I,J,K,:) = ref(I,J,K,:)+dt*(-Fp/eps^2+(ref(1:Nx,J,K,:)+ref(3:Nx+2,J,K,:)+ref(I,1:Ny,K,:)+ref(I,3:Ny+2,K,:)+ref(I,J,1:Nz,:)+ref(I,J,3:Nz+2,:)-6*ref(I,J,K,:))/h^2+beta);
    ref(1,:,:,:) = ref(2,:,:,:); ref(Nx+2,:,:,:) = ref(Nx+1,:,:,:);
    ref(:,1,:,:) = ref(:,2,:,:); ref(:,Ny+2,:,:) = ref(:,Ny+1,:,:);
    ref(:,:,1,:) = ref(:,:,2,:); ref(:,:,Nz+2,:) = ref(:,:,Nz+1,:);
end

L = 1:6;
err = zeros(1,numel(L)); dts = err;
for il = 1:numel(L)
dt = dtmax/2^L(il); dts(il) = dt;
phi = phi0;
for it = 1:16*2^L(il)
    Fp = phi(I,J,K,:).*(phi(I,J,K,:)-0.5).*(phi(I,J,K,:)-1);
    beta = sum(Fp,4)/(N*eps^2);
    phi(I,J,K,:) = phi(I,J,K,:)+dt*(-Fp/eps^2+(phi(1:Nx,J,K,:)+phi(3:Nx+2,J,K,:)+phi(I,1:Ny,K,:)+phi(I,3:Ny+2,K,:)+phi(I,J,1:Nz,:)+phi(I,J,3:Nz+2,:)-6*phi(I,J,K,:))/h^2+beta);
    phi(1,:,:,:) = phi(2,:,:,:); phi(Nx+2,:,:,:) = phi(Nx+1,:,:,:);
    phi(:,1,:,:) = phi(:,2,:,:); phi(:,Ny+2,:,:) = phi(:,Ny+1,:,:);
    phi(:,:,1,:) = phi(:,:,2,:); phi(:,:,Nz+2,:) = phi(:,:,Nz+1,:);
end
err(il) = sqrt(sum((ref(I,J,K,:)-phi(I,J,K,:)).^2,'all')/(N*Nx*Ny*Nz));
end
rate = log2(err(1:end-1)./err(2:end));
fprintf('3D err  = %s\n',sprintf('%.4e ',err));
fprintf('3D rate = %s\n',sprintf('%.4f ',rate));

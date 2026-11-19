clc; clear all; close all;
format long e

Lx = 0; Rx = 1; Ly = 0; Ry = 1; Lz = 0; Rz = 1;
N = 4;
eps = 0.05;
T = 1.0e-3;
MM = [32 64 128 256 512 1024];

err = zeros(1,length(MM)-1);
for k = 1:length(MM)
    Nx = MM(k); h = (Rx-Lx)/Nx; Ny = (Ry-Ly)/h; Nz = (Rz-Lz)/h;
    x = linspace(Lx-0.5*h,Rx+0.5*h,Nx+2);
    y = linspace(Ly-0.5*h,Ry+0.5*h,Ny+2);
    z = linspace(Lz-0.5*h,Rz+0.5*h,Nz+2);
    Nt = 20*(Nx/32)^2; dt = T/Nt;
    [Nx dt 2*eps^2*h^2/(h^2+12*eps^2)]

    [xx,yy,zz] = ndgrid(x,y,z);
    c1 = cos(pi*xx).*cos(pi*yy).*cos(pi*zz); c2 = cos(2*pi*xx).*cos(2*pi*yy).*cos(2*pi*zz);
    phi = cat(4,1/4+0.1*c1,1/4+0.1*c2,1/4-0.1*c1,1/4-0.1*c2);
    phi(1,:,:,:) = phi(2,:,:,:); phi(Nx+2,:,:,:) = phi(Nx+1,:,:,:);
    phi(:,1,:,:) = phi(:,2,:,:); phi(:,Ny+2,:,:) = phi(:,Ny+1,:,:);
    phi(:,:,1,:) = phi(:,:,2,:); phi(:,:,Nz+2,:) = phi(:,:,Nz+1,:);
    for it = 1:Nt
        Fp = phi(2:Nx+1,2:Ny+1,2:Nz+1,:).*(phi(2:Nx+1,2:Ny+1,2:Nz+1,:)-0.5).*(phi(2:Nx+1,2:Ny+1,2:Nz+1,:)-1);
        beta = sum(Fp,4)/(N*eps^2);
        phi(2:Nx+1,2:Ny+1,2:Nz+1,:) = phi(2:Nx+1,2:Ny+1,2:Nz+1,:)+dt*(-Fp/eps^2+(phi(1:Nx,2:Ny+1,2:Nz+1,:)+phi(3:Nx+2,2:Ny+1,2:Nz+1,:) ...
            +phi(2:Nx+1,1:Ny,2:Nz+1,:)+phi(2:Nx+1,3:Ny+2,2:Nz+1,:)+phi(2:Nx+1,2:Ny+1,1:Nz,:)+phi(2:Nx+1,2:Ny+1,3:Nz+2,:)-6*phi(2:Nx+1,2:Ny+1,2:Nz+1,:))/h^2+beta);
        phi(1,:,:,:) = phi(2,:,:,:); phi(Nx+2,:,:,:) = phi(Nx+1,:,:,:);
        phi(:,1,:,:) = phi(:,2,:,:); phi(:,Ny+2,:,:) = phi(:,Ny+1,:,:);
        phi(:,:,1,:) = phi(:,:,2,:); phi(:,:,Nz+2,:) = phi(:,:,Nz+1,:);
    end
    if k > 1
        pf = phi(2:Nx+1,2:Ny+1,2:Nz+1,:);
        pc = (pf(1:2:end,1:2:end,1:2:end,:)+pf(2:2:end,1:2:end,1:2:end,:)+pf(1:2:end,2:2:end,1:2:end,:)+pf(2:2:end,2:2:end,1:2:end,:) ...
            +pf(1:2:end,1:2:end,2:2:end,:)+pf(2:2:end,1:2:end,2:2:end,:)+pf(1:2:end,2:2:end,2:2:end,:)+pf(2:2:end,2:2:end,2:2:end,:))/8;
        err(k-1) = sqrt((2*h)^3*sum((prev-pc).^2,'all')/N);
    end
    prev = phi(2:Nx+1,2:Ny+1,2:Nz+1,:);
end
rate = [NaN log2(err(1:end-1)./err(2:end))];

[MM(1:end-1); err; rate]

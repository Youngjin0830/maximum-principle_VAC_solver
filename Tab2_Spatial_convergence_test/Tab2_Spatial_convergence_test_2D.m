clc; clear all; close all;
format long e

Lx = 0; Rx = 1; Ly = 0; Ry = 1;
N = 4;
eps = 0.05;
T = 1.0e-3;
MM = [32 64 128 256 512 1024];

sol = cell(1,length(MM));
for k = 1:length(MM)
    Nx = MM(k); h = (Rx-Lx)/Nx; Ny = (Ry-Ly)/h;
    x = linspace(Lx-0.5*h,Rx+0.5*h,Nx+2);
    y = linspace(Ly-0.5*h,Ry+0.5*h,Ny+2);
    Nt = 20*(Nx/32)^2; dt = T/Nt;
    [Nx dt 2*eps^2*h^2/(h^2+8*eps^2)]

    [xx,yy] = ndgrid(x,y);
    c1 = cos(pi*xx).*cos(pi*yy); c2 = cos(2*pi*xx).*cos(2*pi*yy);
    phi = cat(3,1/4+0.1*c1,1/4+0.1*c2,1/4-0.1*c1,1/4-0.1*c2);
    phi(1,:,:) = phi(2,:,:); phi(Nx+2,:,:) = phi(Nx+1,:,:);
    phi(:,1,:) = phi(:,2,:); phi(:,Ny+2,:) = phi(:,Ny+1,:);
    for it = 1:Nt
        Fp = phi(2:Nx+1,2:Ny+1,:).*(phi(2:Nx+1,2:Ny+1,:)-0.5).*(phi(2:Nx+1,2:Ny+1,:)-1);
        beta = sum(Fp,3)/(N*eps^2);
        phi(2:Nx+1,2:Ny+1,:) = phi(2:Nx+1,2:Ny+1,:)+dt*(-Fp/eps^2+(phi(1:Nx,2:Ny+1,:)+phi(3:Nx+2,2:Ny+1,:)+phi(2:Nx+1,1:Ny,:)+phi(2:Nx+1,3:Ny+2,:)-4*phi(2:Nx+1,2:Ny+1,:))/h^2+beta);
        phi(1,:,:) = phi(2,:,:); phi(Nx+2,:,:) = phi(Nx+1,:,:);
        phi(:,1,:) = phi(:,2,:); phi(:,Ny+2,:) = phi(:,Ny+1,:);
    end
    sol{k} = phi(2:Nx+1,2:Ny+1,:);
end

err = zeros(1,length(MM)-1);
for k = 1:length(MM)-1
    h = (Rx-Lx)/MM(k);
    pf = sol{k+1};
    pc = (pf(1:2:end,1:2:end,:)+pf(2:2:end,1:2:end,:)+pf(1:2:end,2:2:end,:)+pf(2:2:end,2:2:end,:))/4;
    err(k) = sqrt(h^2*sum((sol{k}-pc).^2,'all')/N);
end
rate = [NaN log2(err(1:end-1)./err(2:end))];

[MM(1:end-1); err; rate]

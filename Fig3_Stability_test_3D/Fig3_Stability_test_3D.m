clc; clear all; close all;
fig = figure(1);

Nx = 128;
Lx = 0; Rx = 1; h = (Rx-Lx)/Nx;
Ly = 0; Ry = 1; Ny = (Ry-Ly)/h;
Lz = 0; Rz = 1; Nz = (Rz-Lz)/h;
x = linspace(Lx-0.5*h,Rx+0.5*h,Nx+2);
y = linspace(Ly-0.5*h,Ry+0.5*h,Ny+2);
z = linspace(Lz-0.5*h,Rz+0.5*h,Nz+2);
N = 4;
[xi,yi,zi] = meshgrid(x(2:Nx+1),y(2:Ny+1),z(2:Nz+1));

m = 4;
eps = m*h/(4*sqrt(2)*atanh(0.9));

dtmax = 2*eps^2*h^2/(h^2+12*eps^2);

rng(33)
w = 1+0.1*(2*rand(Nx,Ny,Nz,N)-1);
phi0 = zeros(Nx+2,Ny+2,Nz+2,N);
phi0(2:Nx+1,2:Ny+1,2:Nz+1,:) = w./sum(w,4);
clear w
phi0(1,:,:,:) = phi0(2,:,:,:); phi0(Nx+2,:,:,:) = phi0(Nx+1,:,:,:);
phi0(:,1,:,:) = phi0(:,2,:,:); phi0(:,Ny+2,:,:) = phi0(:,Ny+1,:,:);
phi0(:,:,1,:) = phi0(:,:,2,:); phi0(:,:,Nz+2,:) = phi0(:,:,Nz+1,:);

cl = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
fs = 21; ms = 10;
I = 2:Nx+1; J = 2:Ny+1; K = 2:Nz+1;

fac = [1 1.25];
Ntl = [1920 1536];
itsnap = {[240 600 1920], [192 480 1536]};
tags = {'a','b'};
for c = 1:2
dt = fac(c)*dtmax;
Nt = Ntl(c);
phi = phi0;

maxi = zeros(Nt+1,1);
mini = maxi;
ene = maxi;
sume = maxi;
maxi(1) = max(phi(I,J,K,:),[],'all');
mini(1) = min(phi(I,J,K,:),[],'all');
sume(1) = max(abs(sum(phi(I,J,K,:),4)-1),[],'all');
ene(1) = sum(h^3*phi(I,J,K,:).^2.*(1-phi(I,J,K,:)).^2/(4*eps^2),'all')+h*(sum((phi(3:Nx+1,J,K,:)-phi(2:Nx,J,K,:)).^2,'all')+sum((phi(I,3:Ny+1,K,:)-phi(I,2:Ny,K,:)).^2,'all')+sum((phi(I,J,3:Nz+1,:)-phi(I,J,2:Nz,:)).^2,'all'))/2;

js = 1;
for it = 1:Nt

Fp = phi(I,J,K,:).*(phi(I,J,K,:)-0.5).*(phi(I,J,K,:)-1);
beta = sum(Fp,4)/(N*eps^2);
phi(I,J,K,:) = phi(I,J,K,:)+dt*(-Fp/eps^2+(phi(1:Nx,J,K,:)+phi(3:Nx+2,J,K,:)+phi(I,1:Ny,K,:)+phi(I,3:Ny+2,K,:)+phi(I,J,1:Nz,:)+phi(I,J,3:Nz+2,:)-6*phi(I,J,K,:))/h^2+beta);

phi(1,:,:,:) = phi(2,:,:,:); phi(Nx+2,:,:,:) = phi(Nx+1,:,:,:);
phi(:,1,:,:) = phi(:,2,:,:); phi(:,Ny+2,:,:) = phi(:,Ny+1,:,:);
phi(:,:,1,:) = phi(:,:,2,:); phi(:,:,Nz+2,:) = phi(:,:,Nz+1,:);

maxi(it+1) = max(phi(I,J,K,:),[],'all');
mini(it+1) = min(phi(I,J,K,:),[],'all');
sume(it+1) = max(abs(sum(phi(I,J,K,:),4)-1),[],'all');
ene(it+1) = sum(h^3*phi(I,J,K,:).^2.*(1-phi(I,J,K,:)).^2/(4*eps^2),'all')+h*(sum((phi(3:Nx+1,J,K,:)-phi(2:Nx,J,K,:)).^2,'all')+sum((phi(I,3:Ny+1,K,:)-phi(I,2:Ny,K,:)).^2,'all')+sum((phi(I,J,3:Nz+1,:)-phi(I,J,2:Nz,:)).^2,'all'))/2;

if js <= 3 && it == itsnap{c}(js)
    clf(fig); hold on; box on; grid on;
    set(gca,'fontsize',fs);
    set(gca, 'TickLabelInterpreter', 'latex');
    set(gcf,'PaperPositionMode','auto','Position',[500 500 450 400])
    for p = 1:N
        V = permute(phi(I,J,K,p),[2 1 3]);
        patch(isosurface(xi,yi,zi,V,0.5),'facecolor',cl(p,:),'facealpha',0.6,'edgecolor','none');
        patch(isocaps(xi,yi,zi,V,0.5,'above'),'facecolor',cl(p,:),'facealpha',0.6,'edgecolor','none');
    end
    axis image
    view([-30 15]);
    camlight('right'); daspect([1 1 1]); lighting phong;
    axis([Lx Rx Ly Ry Lz Rz])
    text('Interpreter','latex','String','$x$','Position',[0.704 -0.267 -0.016],'FontSize',fs)
    text('Interpreter','latex','String','$y$','Position',[-0.148 0.796 -0.021],'FontSize',fs-2)
    text('Interpreter','latex','String','$z$','Position',[-0.274 0.791 0.868],'FontSize',fs)
    print(fig,'-depsc',sprintf('figb3_%s%d.eps',tags{c},js))
    js = js+1;
end

end

fprintf('3D dt = %.2f dtmax: steps = %d, Max = %.6e, Min = %.6e, max|sum-1| = %.3e, energy increases = %d\n', ...
    fac(c), Nt, max(maxi), min(mini), max(sume), sum(diff(ene) > 1e-12*ene(1)));
t{c} = (0:Nt)'*dt; mx{c} = maxi; mn{c} = mini; E{c} = ene;
end

n1 = Ntl(1)+1; n2 = Ntl(2)+1;
clf(fig); hold on; box on; grid on;
set(gca,'fontsize',fs);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto','Position',[100 100 780 400])
plot(t{1},mx{1},'b-^','markersize',ms,'color','#6FA8DC','markerfacecolor','none','linewidth',1.5,'MarkerIndices',1:200:n1)
plot(t{1},mn{1},'b-v','markersize',ms,'color','#6FA8DC','markerfacecolor','none','linewidth',1.5,'MarkerIndices',1:200:n1)
plot(t{2},mx{2},'r:^','markersize',ms,'color','#E77C73','markerfacecolor','#E77C73','linewidth',2,'MarkerIndices',100:200:n2)
plot(t{2},mn{2},'r:v','markersize',ms,'color','#E77C73','markerfacecolor','#E77C73','linewidth',2,'MarkerIndices',100:200:n2)
axis([0 t{1}(end) -0.5 1.5])
text('Interpreter','latex','String','$t$','FontSize',24,'Units','normalized','Position',[0.975 -0.073])
leg = legend('$Max(\mbox{\boldmath$\phi$}^n),~\Delta t=\Delta t_{max}$','$Min(\mbox{\boldmath$\phi$}^n),~\Delta t=\Delta t_{max}$','$Max(\mbox{\boldmath$\phi$}^n),~\Delta t=1.25\Delta t_{max}$','$Min(\mbox{\boldmath$\phi$}^n),~\Delta t=1.25\Delta t_{max}$');
set(leg,'Interpreter','latex','Location','east','NumColumns',2,'FontSize',15)
exportgraphics(gcf,'figb3_c.eps','ContentType','vector','Padding',3)

e1 = E{1}/E{1}(2); e2 = E{2}/E{2}(2);
emax = max([e1;e2]);
elo = 10^floor(log10(min([e1;e2])));
ek = log10(elo):floor(log10(emax));
clf(fig);
loglog(t{1},e1,'b-','linewidth',1.5,'color','#6FA8DC'); hold on; box on; grid on;
set(gca,'fontsize',fs);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto','Position',[100 100 740 400])
loglog(t{2},e2,'r:','color','#E77C73','linewidth',2)
axis([t{1}(2) t{1}(end) elo emax])
yticks(10.^ek); yticklabels(compose('$10^{%d}$',ek));
xticks(10.^[-5 -4 -3 -2])
text('Interpreter','latex','String','$t$','FontSize',24,'Units','normalized','Position',[0.795714285714286 -0.0784804928131417])
text('Interpreter','latex','String','$\frac{\mathcal{E}_h(\mbox{\boldmath$\phi$}^n)}{\mathcal{E}_h(\mbox{\boldmath$\phi$}^1)}$','FontSize',20,'Units','normalized','Position',[-0.153554006968641 0.875954825462013])
leg = legend('$\Delta t=\Delta t_{max}$','$\Delta t=1.25\Delta t_{max}$');
set(leg,'Interpreter','latex','Location','southwest')
print(fig,'-depsc','figb3_d.eps')

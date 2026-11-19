clc; clear all; close all;
fig = figure(1);

Nx = 64;
Lx = 0; Rx = 1; h = (Rx-Lx)/Nx;
Ly = 0; Ry = 1; Ny = (Ry-Ly)/h;

m = 4;
eps = m*h/(4*sqrt(2)*atanh(0.9));

dtmax = 2*eps^2*h^2/(h^2+8*eps^2);
dt = dtmax;
Nt = 100;
T = Nt*dt;

Nlist = 2:7;
maxi = zeros(Nt+1,numel(Nlist));
mini = maxi;
sume = maxi;
for jN = 1:numel(Nlist)
N = Nlist(jN);
rng(10+N)
w = 0.1*rand(Nx,Ny,N);
q = randi(N,Nx,Ny);
for p = 1:N
    wp = w(:,:,p); wp(q == p) = 1; w(:,:,p) = wp;
end
phi = zeros(Nx+2,Ny+2,N);
phi(2:Nx+1,2:Ny+1,:) = w./sum(w,3);
phi(1,:,:) = phi(2,:,:); phi(Nx+2,:,:) = phi(Nx+1,:,:);
phi(:,1,:) = phi(:,2,:); phi(:,Ny+2,:) = phi(:,Ny+1,:);
maxi(1,jN) = max(phi(2:Nx+1,2:Ny+1,:),[],'all');
mini(1,jN) = min(phi(2:Nx+1,2:Ny+1,:),[],'all');
sume(1,jN) = max(abs(sum(phi(2:Nx+1,2:Ny+1,:),3)-1),[],'all');

for it = 1:Nt
Fp = phi(2:Nx+1,2:Ny+1,:).*(phi(2:Nx+1,2:Ny+1,:)-0.5).*(phi(2:Nx+1,2:Ny+1,:)-1);
beta = sum(Fp,3)/(N*eps^2);
phi(2:Nx+1,2:Ny+1,:) = phi(2:Nx+1,2:Ny+1,:)+dt*(-Fp/eps^2+(phi(1:Nx,2:Ny+1,:)+phi(3:Nx+2,2:Ny+1,:)+phi(2:Nx+1,1:Ny,:)+phi(2:Nx+1,3:Ny+2,:)-4*phi(2:Nx+1,2:Ny+1,:))/h^2+beta);
phi(1,:,:) = phi(2,:,:); phi(Nx+2,:,:) = phi(Nx+1,:,:);
phi(:,1,:) = phi(:,2,:); phi(:,Ny+2,:) = phi(:,Ny+1,:);
maxi(it+1,jN) = max(phi(2:Nx+1,2:Ny+1,:),[],'all');
mini(it+1,jN) = min(phi(2:Nx+1,2:Ny+1,:),[],'all');
sume(it+1,jN) = max(abs(sum(phi(2:Nx+1,2:Ny+1,:),3)-1),[],'all');
end
fprintf('N = %d: Max = %.15f, Min = %.6e, Max(T) = %.4f, Min(T) = %.4e, max|sum-1| = %.3e\n',N,max(maxi(:,jN)),min(mini(:,jN)),maxi(end,jN),mini(end,jN),max(sume(:,jN)));
end

cl = lines(6);
mk = {'o','s','^','v','d','>'};
fs = 11; ms = 5;
mi = unique(round(logspace(0,log10(Nt),10)))+1;
tt = (0:Nt)'*dt;
clf(fig); hold on; box on; grid on;
set(gca,'fontsize',fs);
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf,'PaperPositionMode','auto','Position',[100 100 780 270])
for jN = 1:numel(Nlist)
    plot(tt,maxi(:,jN),'-','color',cl(jN,:),'marker',mk{jN},'markersize',ms,'linewidth',1.2,'MarkerIndices',mi)
end
for jN = 1:numel(Nlist)
    plot(tt,mini(:,jN),'-','color',cl(jN,:),'marker',mk{jN},'markersize',ms,'linewidth',1.2,'MarkerIndices',mi,'HandleVisibility','off')
end
set(gca,'XScale','log')
axis([tt(2) T -0.1 1.1])
text('Interpreter','latex','String','$t$','FontSize',13,'Units','normalized','Position',[1.01 -0.06])
leg = legend('$N=2$','$N=3$','$N=4$','$N=5$','$N=6$','$N=7$');
set(leg,'Interpreter','latex','Location','northoutside','Orientation','horizontal','NumColumns',3,'FontSize',10)
lp = get(leg,'Position'); ap = get(gca,'Position');
atop = 0.97-lp(4);
set(gca,'Position',[ap(1) ap(2) ap(3) atop-ap(2)])
set(leg,'Position',[lp(1) atop+0.01 lp(3) lp(4)])
exportgraphics(gcf,'fig3.eps','ContentType','vector','Padding',3)

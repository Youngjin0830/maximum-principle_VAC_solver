clc; clear all; close all;
fig = figure(1);

Nx = 128;
Lx = 0; Rx = 1; h = (Rx-Lx)/Nx;
x = linspace(Lx-0.5*h,Rx+0.5*h,Nx+2);
N = 4;

m = 4;
eps = m*h/(4*sqrt(2)*atanh(0.9));

dtmax = 2*eps^2*h^2/(h^2+4*eps^2);

rng(31)
w = 1+0.1*(2*rand(Nx,N)-1);
phi0 = zeros(Nx+2,N);
phi0(2:Nx+1,:) = w./sum(w,2);
phi0(1,:) = phi0(2,:); phi0(Nx+2,:) = phi0(Nx+1,:);

cl = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
mk = {'o','s','^','d'};
fs = 21; ms = 10;

fac = [1 1.5];
Ntl = [1920 1280];
itsnap = {[120 600 1920], [80 400 1280]};
tags = {'a','b'};
for c = 1:2
dt = fac(c)*dtmax;
Nt = Ntl(c);
phi = phi0;

maxi = zeros(Nt+1,1);
mini = maxi;
ene = maxi;
sume = maxi;
maxi(1) = max(phi(2:Nx+1,:),[],'all');
mini(1) = min(phi(2:Nx+1,:),[],'all');
sume(1) = max(abs(sum(phi(2:Nx+1,:),2)-1),[],'all');
ene(1) = sum(h*phi(2:Nx+1,:).^2.*(1-phi(2:Nx+1,:)).^2/(4*eps^2),'all')+sum((phi(3:Nx+1,:)-phi(2:Nx,:)).^2,'all')/(2*h);

js = 1;
for it = 1:Nt

Fp = phi(2:Nx+1,:).*(phi(2:Nx+1,:)-0.5).*(phi(2:Nx+1,:)-1);
beta = sum(Fp,2)/(N*eps^2);
phi(2:Nx+1,:) = phi(2:Nx+1,:)+dt*(-Fp/eps^2+(phi(1:Nx,:)+phi(3:Nx+2,:)-2*phi(2:Nx+1,:))/h^2+beta);

phi(1,:) = phi(2,:); phi(Nx+2,:) = phi(Nx+1,:);

maxi(it+1) = max(phi(2:Nx+1,:),[],'all');
mini(it+1) = min(phi(2:Nx+1,:),[],'all');
sume(it+1) = max(abs(sum(phi(2:Nx+1,:),2)-1),[],'all');
ene(it+1) = sum(h*phi(2:Nx+1,:).^2.*(1-phi(2:Nx+1,:)).^2/(4*eps^2),'all')+sum((phi(3:Nx+1,:)-phi(2:Nx,:)).^2,'all')/(2*h);

if js <= 3 && it == itsnap{c}(js)
    clf(fig); hold on; box on; grid on;
    set(gca,'fontsize',fs);
    set(gca, 'TickLabelInterpreter', 'latex');
    set(gcf,'PaperPositionMode','auto','Position',[500 500 450 400])
    for p = 1:N
        plot(x(2:Nx+1),phi(2:Nx+1,p),'-','color',cl(p,:),'marker',mk{p},'markersize',7,'MarkerIndices',3*p-2:12:Nx,'linewidth',1.5)
    end
    axis([Lx Rx -0.5 1.5])
    text('Interpreter','latex','String','$x$','FontSize',24,'Units','normalized','Position',[0.735 -0.058])
    text('Interpreter','latex','String','$\mbox{\boldmath$\phi$}$','FontSize',24,'Units','normalized','Position',[-0.096 0.883])
    if c == 1 && js == 1
        legend('$\phi_1$','$\phi_2$','$\phi_3$','$\phi_4$','Interpreter','latex','Location','north','NumColumns',4,'FontSize',14)
    end
    print(fig,'-depsc',sprintf('figb1_%s%d.eps',tags{c},js))
    js = js+1;
end

end

fprintf('1D dt = %.2f dtmax: steps = %d, Max = %.6e, Min = %.6e, max|sum-1| = %.3e, energy increases = %d\n', ...
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
leg = legend('$Max(\mbox{\boldmath$\phi$}^n),~\Delta t=\Delta t_{max}$','$Min(\mbox{\boldmath$\phi$}^n),~\Delta t=\Delta t_{max}$','$Max(\mbox{\boldmath$\phi$}^n),~\Delta t=1.5\Delta t_{max}$','$Min(\mbox{\boldmath$\phi$}^n),~\Delta t=1.5\Delta t_{max}$');
set(leg,'Interpreter','latex','Location','east','NumColumns',2,'FontSize',15)
exportgraphics(gcf,'figb1_c.eps','ContentType','vector','Padding',3)

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
text('Interpreter','latex','String','$t$','FontSize',24,'Units','normalized','Position',[0.976898954703833 -0.0846406570841889])
text('Interpreter','latex','String','$\frac{\mathcal{E}_h(\mbox{\boldmath$\phi$}^n)}{\mathcal{E}_h(\mbox{\boldmath$\phi$}^1)}$','FontSize',20,'Units','normalized','Position',[-0.15006968641115 0.95911704312115])
leg = legend('$\Delta t=\Delta t_{max}$','$\Delta t=1.5\Delta t_{max}$');
set(leg,'Interpreter','latex','Location','southwest')
print(fig,'-depsc','figb1_d.eps')

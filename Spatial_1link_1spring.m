

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           Parameter Definitions                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k1      = rand;
m1      = rand;
g       = 9.81;

xp1 = rand;
xs1 = rand;

yp1 = rand;
ys1 = 0;

zp1 = rand;
zs1 = 0;

xc1  = (k1*xp1*xs1)/(g*m1);
yc1  = (k1*yp1*xs1)/(g*m1);
zc1  = (k1*zp1*xs1)/(g*m1);

Eq =    [   k1*xp1*xs1 - g*m1*xc1
            k1*xp1*ys1
            k1*xp1*zs1
            k1*yp1*xs1 - g*m1*yc1
            k1*yp1*ys1
            k1*yp1*zs1
            k1*zp1*xs1 - g*m1*zc1
            k1*zp1*ys1
            k1*zp1*zs1];



Rz = @(a) [ cos(a) -sin(a)  0;
            sin(a)  cos(a)  0;
            0       0       1];

Ry = @(a) [ cos(a)  0       sin(a);
            0       1       0;
            -sin(a) 0       cos(a)];

Rx = @(a) [ 1       0       0 ;
            0       cos(a)  -sin(a);
            0       sin(a)  cos(a);];

Rzyx = @(a,b,c) Rz(a)*Ry(b)*Rx(c);

H = @(a,b,c,t) [Rzyx(a,b,c), [t;0;0] ; [0,0,0,1]];

P = [xp1;yp1;zp1;1];
S = [xs1;ys1;zs1;1];
C = [xc1;yc1;zc1;1];

alpha1 = 0:0.01:2*pi;
alpha2 = sqrt(2).*alpha1;
alpha3 = pi.*alpha1;

for k = 1:numel(alpha1)
    P_prime(:,k) = H(alpha1(k),alpha2(k),alpha3(k),0)*P;
    S_prime(:,k) = S;
    M_prime(:,k) = H(alpha1(k),alpha2(k),alpha3(k),0)*C;
end

D = P_prime-S_prime;

for k = 1:numel(alpha1)
    
    Ue(k) = 0.5*k1*(D(:,k).'*D(:,k));
    Ug(k) = m1*g*M_prime(1,k);
    
end

Ug = Ug - min(Ug);

Utotal = Ue+Ug;

Fsize = 16;


% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');

figure('color',[1,1,1])
set(gca,'TickLabelInterpreter', 'latex','fontsize',Fsize);
hold on
xlabel('$\alpha$ [rad]','interpreter','latex','fontsize',Fsize)
xticks([0,0.5*pi,1*pi,1.5*pi,2*pi])
xticklabels({'$0$','$0.5\pi$','$\pi$','$1.5\pi$','$2\pi$'})
ylabel('Energy [Nm]','interpreter','latex','fontsize',Fsize)
grid on


plot(alpha1,Utotal,'-k','linewidth',2);hold on
plot(alpha1,Ug,'--k','linewidth',2)
plot(alpha1,Ue,'-r','linewidth',2)
plot(2*pi,0.01,'w.')


legend({'Total $U$','Gravitational $U_g$', 'Spring $U_e$'},'interpreter','latex','location','Northwest','fontsize',Fsize)











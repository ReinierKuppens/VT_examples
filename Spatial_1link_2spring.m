




clear all 
close all 
clc 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           Parameter Definitions                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k1      = rand;
k2      = rand; 
m1      = rand;
g       = 9.81;


zp1 = rand; 
xs1 = rand;
xp2 = rand;
yp2 = rand;
zp2 = rand;
xs2 = rand;
ys2 = rand;
zs2 = rand;

xc1     = xp2*(k1*zp1*xs1+k2*zp2*xs2)/(g*m1*zp2); 
yc1     = yp2*(k1*zp1*xs1+k2*zp2*xs2)/(g*m1*zp2);
zc1     = zp2*(k1*zp1*xs1+k2*zp2*xs2)/(g*m1*zp2);

xp1     = (xp2*zp1)/zp2;
yp1     = (yp2*zp1)/zp2;

ys1     = -(k2*zp2*ys2)/(k1*zp1);
zs1     = -(k2*zp2*zs2)/(k1*zp1);

 P1 = [xp1;yp1;zp1;1];
 P2 = [xp2;yp2;zp2;1];
 S1 = [xs1;ys1;zs1;1];
 S2 = [xs2;ys2;zs2;1];

 C = [xc1;yc1;zc1;1];


Eq =    [    k1*xp1*xs1 + k2*xp2*xs2 - g*m1*xc1
             k1*xp1*ys1 + k2*xp2*ys2
             k1*xp1*zs1 + k2*xp2*zs2
             k1*yp1*xs1 + k2*yp2*xs2 - g*m1*yc1
             k1*yp1*ys1 + k2*yp2*ys2
             k1*yp1*zs1 + k2*yp2*zs2
             k1*zp1*xs1 + k2*zp2*xs2 - g*m1*zc1
             k1*zp1*ys1 + k2*zp2*ys2
             k1*zp1*zs1 + k2*zp2*zs2];


k1      = rand;
k2      = rand; 

m1      = rand;
g       = 9.81;

P1 = rand([3,1])
P2 = rand([3,1])

xs1 = rand; 
xs2 = rand; 

C = (k1*xs1/(g*m1))*P1 + (k2*xs2/(g*m1))*P2

zs1 = rand;
ys2 = rand;
zs2 = rand; 


Zs1Ys2 = zs1*ys2;
ys1 = Zs1Ys2/zs2;
         
zs1 = 0;
ys1 = 0;
zs2 = 0;
ys2 = 0;
   
P1 = [P1;1];
P2 = [P2;1];
S1 = [xs1;ys1;zs1;1];
S2 = [xs2;ys2;zs2;1];
C  = [C;1];






   Rz = @(a) [cos(a) -sin(a)  0; 
            sin(a)  cos(a)  0; 
            0       0       1];
        
Ry = @(a) [cos(a)  0       sin(a); 
            0       1       0; 
            -sin(a) 0       cos(a)];
        
 Rx = @(a) [1       0       0 ; 
            0       cos(a)  -sin(a); 
            0       sin(a)  cos(a);];
 Rzyx = @(a,b,c) Rz(a)*Ry(b)*Rx(c);
 
 H = @(a,b,c,t) [Rzyx(a,b,c), [t;0;0] ; [0,0,0,1]];
 alpha1 = 0:0.01:2*pi;
 alpha2 = sqrt(2).*alpha1;
 alpha3 = pi.*alpha1;
 
 
 for k = 1:numel(alpha1) 
%     keyboard 
     P1_prime(:,k) = H(alpha1(k),alpha2(k),alpha3(k),0)*P1;
     S1_prime(:,k) = S1;
     P2_prime(:,k) = H(alpha1(k),alpha2(k),alpha3(k),0)*P2;
     S2_prime(:,k) = S2;
     M_prime(:,k) = H(alpha1(k),alpha2(k),alpha3(k),0)*C;
          
 end
 
 D1 = P1_prime-S1_prime;
 D2 = P2_prime-S2_prime; 
 
 for k = 1:numel(alpha1) 
    
     Ue1(k) = 0.5*k1*(D1(:,k).'*D1(:,k));
     Ue2(k) = 0.5*k2*(D2(:,k).'*D2(:,k));

     Ug(k) = m1*g*M_prime(1,k);
     
 end
 
 Ug = Ug - min(Ug);
 
 Utotal = Ue1+Ue2+Ug;
 
 %%
 
 Fsize = 16; 
figure('color',[1,1,1])
set(gca,'TickLabelInterpreter', 'latex','fontsize',Fsize);
hold on
xlabel('$\alpha$ [rad]','interpreter','latex','fontsize',Fsize)
xticks([0,0.5*pi,1*pi,1.5*pi,2*pi])
xticklabels({'$0$','$0.5\pi$','$\pi$','$1.5\pi$','$2\pi$'})
ylabel('Energy [Nm]','interpreter','latex','fontsize',Fsize)
grid on

plot(alpha1,Utotal,'-k','linewidth',2)
 plot(alpha1,Ug,'--k','linewidth',2)

 plot(alpha1,Ue1,'-r','linewidth',2)
 plot(alpha1,Ue2,'-r','linewidth',2)
 plot(2*pi,0.01,'w.')
     
legend({'Total $U$','Gravitational $U_g$', 'Spring $U_e$'},'interpreter','latex','location','Northwest','fontsize',Fsize) 




x = -10:0.01:10;E1 = P1+x.*(P2-P1);

 E2 = x.*1.*(C);

% figure;plot3(E1(1,:),E1(2,:),E1(3,:),'k');hold on; 
% plot3(E2(1,:),E2(2,:),E2(3,:),'r')

%%
% figure;hold on 
% plot3([xp1;xs1],[yp1;ys1],[zp1;zs1],'linewidth',2,'color','r')
% plot3([xp2;xs2],[yp2;ys2],[zp2;zs2],'linewidth',2,'color','b')
% plot3([0;xc1],[0;yc1],[0;zc1],'linewidth',2,'color','k')
% plot3(0,0,0,'marker','.','markersize',20)
% 
% plot3([xs1,xs2],[ys1,ys2],[zs1,zs2],'linewidth',2,'color','y')
% plot3([0,2],[0,0],[0,0],'linewidth',2,'color','c')




clear all 
close all 
clc 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           Parameter Definitions                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSpring = 5; 

m = 0.5*rand;
g = 9.81;

kn2 = rand([1,nSpring-2]);

xpn2 = rand([1,nSpring-2]);
ypn2 = rand([1,nSpring-2]);
zpn2 = rand([1,nSpring-2]);

xsn2 = rand([1,nSpring-2]);
ysn2 = rand([1,nSpring-2]);
zsn2 = rand([1,nSpring-2]);

k1  = rand;
xs1 = rand;
ys1 = rand;
zs1 = rand;

k2  = rand;
xs2 = rand;
ys2 = rand;
zs2 = rand;

C1 = sum(kn2.*xpn2.*xsn2);
C2 = sum(kn2.*xpn2.*ysn2);
C3 = sum(kn2.*xpn2.*zsn2);
C4 = sum(kn2.*ypn2.*xsn2);
C5 = sum(kn2.*ypn2.*ysn2);
C6 = sum(kn2.*ypn2.*zsn2);
C7 = sum(kn2.*zpn2.*xsn2);
C8 = sum(kn2.*zpn2.*ysn2);
C9 = sum(kn2.*zpn2.*zsn2);



xc1 = (C1*ys1*zs2-C1*ys2*zs1-C2*xs1*zs2+C2*xs2*zs1+C3*xs1*ys2-C3*xs2*ys1)/(g*m*(ys1*zs2-ys2*zs1));
yc1 = (C4*ys1*zs2-C4*ys2*zs1-C5*xs1*zs2+C5*xs2*zs1+C6*xs1*ys2-C6*xs2*ys1)/(g*m*(ys1*zs2-ys2*zs1));
zc1 = (C7*ys1*zs2-C7*ys2*zs1-C8*xs1*zs2+C8*xs2*zs1+C9*xs1*ys2-C9*xs2*ys1)/(g*m*(ys1*zs2-ys2*zs1));
xp1 = -(C2*zs2-C3*ys2)/((ys1*zs2-ys2*zs1)*k1);
xp2 = (C2*zs1-C3*ys1)/(k2*(ys1*zs2-ys2*zs1));
yp1 = -(C5*zs2-C6*ys2)/((ys1*zs2-ys2*zs1)*k1);
yp2 = (C5*zs1-C6*ys1)/(k2*(ys1*zs2-ys2*zs1));
zp1 = -(C8*zs2-C9*ys2)/((ys1*zs2-ys2*zs1)*k1);
zp2 = (C8*zs1-C9*ys1)/(k2*(ys1*zs2-ys2*zs1));


Eq = [  C1 + k1*xp1*xs1+k2*xp2*xs2-g*m*xc1
        C2 + k1*xp1*ys1+k2*xp2*ys2
        C3 + k1*xp1*zs1+k2*xp2*zs2
        C4 + k1*yp1*xs1+k2*yp2*xs2-g*m*yc1
        C5 + k1*yp1*ys1+k2*yp2*ys2
        C6 + k1*yp1*zs1+k2*yp2*zs2
        C7 + k1*zp1*xs1+k2*zp2*xs2-g*m*zc1
        C8 + k1*zp1*ys1+k2*zp2*ys2
        C9 + k1*zp1*zs1+k2*zp2*zs2]



    xp = [xp1 xp2 xpn2];
    yp = [yp1 yp2 ypn2];
    zp = [zp1 zp2 zpn2];
    xs = [xs1 xs2 xsn2];
    ys = [ys1 ys2 ysn2];
    zs = [zs1 zs2 zsn2];
    
    Kstif  = [k1 k2 kn2];
    
 P = [xp;yp;zp;ones(1,nSpring)];
 S = [xs;ys;zs;ones(1,nSpring)];


 C = [xc1;yc1;zc1;1];
 
 alpha1 = 0:0.01:2*pi;
 alpha2 = sqrt(2).*alpha1;
 alpha3 = pi.*alpha1;
    
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
 
 for k = 1:numel(alpha1) 
     P_prime(:,:,k) = H(alpha1(k),alpha2(k),alpha3(k),0)*P;
     S_prime(:,:,k) = S;
     M_prime(:,k) = H(alpha1(k),alpha2(k),alpha3(k),0)*C;
 end
 
 D = P_prime-S_prime;
 
  for k = 1:numel(alpha1) 
      for ii = 1:nSpring
        Ue(ii,k) = 0.5*Kstif(ii)*(D(:,ii,k).'*D(:,ii,k));
      end
     Ug(1,k) = m*g*M_prime(1,k);
 end
 
 Ug     = Ug - min(Ug);
 Utotal = sum([Ue;Ug]);
 
Fsize = 16;

figure('color',[1,1,1])
set(gca,'TickLabelInterpreter', 'latex','fontsize',Fsize);
hold on
xlabel('$\alpha$ [rad]','interpreter','latex','fontsize',Fsize)
xticks([0,0.5*pi,1*pi,1.5*pi,2*pi])
xticklabels({'$0$','$0.5\pi$','$\pi$','$1.5\pi$','$2\pi$'})
ylabel('Energy [Nm]','interpreter','latex','fontsize',Fsize)
grid on

 
 plot(alpha1,Utotal,'k','linewidth',2) ;hold on 
 plot(alpha1,Ug,'--k','linewidth',2) 

for k = 1:nSpring
   plot(alpha1,Ue(k,:),'-r','linewidth',2) 
end
plot(2*pi,0.01,'w.')
     
legend({'Total $U$','Gravitational $U_g$', 'Spring $U_e$'},'interpreter','latex','location','Northwest','fontsize',Fsize) 













 
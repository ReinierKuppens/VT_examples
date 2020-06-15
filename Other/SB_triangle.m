

clear all 
close all 
clc 


L = 3;


alpha = 0:0.01:2*pi;
gamma =  pi/6;
zeta  =  pi/3; 
beta  =  zeta - alpha; 

x1 = tan(gamma).*sin(alpha).*L + cos(alpha).*L;
x2 = (sin(alpha).*L)./cos(gamma);
x3 = (sin(beta).*L )./sin(zeta);  

p1 = [  x1;
        zeros(size(x1))];
p2 = [  -cos(2*zeta)*x2;
        -sin(2*zeta)*x2];
p3 = [  cos(zeta)*x3;
        sin(zeta)*x3];
    
    


% figure;plot(alpha,x1)

figure;hold on ; axis equal; 
axis([-L*3,L*3,-L*3,L*3])

plot(cos(0).*[-3*L:0.01:3*L],sin(0).*[-3*L:0.01:3*L])
plot(cos(zeta).*[-3*L:0.01:3*L],sin(zeta).*[-3*L:0.01:3*L])
plot(cos(2*zeta).*[-3*L:0.01:3*L],sin(2*zeta).*[-3*L:0.01:3*L])

R2 = (L/2)/tan(pi/3) + cos(pi/6)*L
R1 = R2/2;

c1x = cos(alpha).*R1;
c1y = sin(alpha).*R1;

midleC = (p1+p2+p3)./3;

plot(cos(alpha).*R2,sin(alpha).*R2,'r')
H2=plot(c1x + midleC(1,1),c1y + midleC(2,1));


% keyboard 
% for k = 1:numel(alpha)
    k=1;
H = plot([p1(1,k),p2(1,k),p3(1,k),p1(1,k)],[p1(2,k),p2(2,k),p3(2,k),p1(2,k)],'k');
% pause(0.01)
% end

% keyboard 

for k = 1:numel(alpha)
    set(H,'Xdata',[p1(1,k),p2(1,k),p3(1,k),p1(1,k)],'Ydata',[p1(2,k),p2(2,k),p3(2,k),p1(2,k)])
    set(H2,'Xdata',c1x + midleC(1,k),'Ydata',c1y + midleC(2,k))
    pause(0.01)
end

%%% check 

dp12    = p1-p2;
dp13    = p1-p3;
dp23    = p2-p3;

np12    = (dp12(1,:).^2 + dp12(2,:).^2).^0.5;
np13    = (dp13(1,:).^2 + dp13(2,:).^2).^0.5;
np23    = (dp23(1,:).^2 + dp23(2,:).^2).^0.5;

sum(diff(np12))
sum(diff(np13))
sum(diff(np23))

k1 = 1; 
U  = 0.5*k1*x1.^2 +  0.5*k1*x2.^2 +  0.5*k1*x3.^2 ;

figure;plot(alpha,U)

qss = [midleC.',alpha.'];





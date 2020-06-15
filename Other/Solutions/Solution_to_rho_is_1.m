
clear all 
close all 
clc 


nSpring  = 5;

rout     = 2;
rin      = 1;

rho     =  rout/rin;

alpha   = linspace(0,2*pi,200);
beta    = -rho*linspace(0,2*pi,200);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Sove the equation for coupling: sum krq*[cos(theta-phi);sin(theta-phi)]

%%%%%%%%%%%%%% 1

rbup = rin;

q_n1 =  2*rbup.*rand([1,nSpring-1]);
r_n1 =  2*rbup.*rand([1,nSpring-1]);
k_n1 =  2*rbup.*rand([1,nSpring-1]);

theta_n1 = 2*pi.*rand([1,nSpring-1]);
phi_n1   = 2*pi.*rand([1,nSpring-1]);


C1      = sum(k_n1.*r_n1.*cos(theta_n1)) - sum(k_n1.*q_n1.*cos(phi_n1));
C2      = sum(k_n1.*r_n1.*sin(theta_n1)) + sum(k_n1.*q_n1.*sin(phi_n1));
C3      = sum(2.*k_n1.*r_n1.*q_n1.*cos(theta_n1-phi_n1));
C4      = sum(2.*k_n1.*r_n1.*q_n1.*sin(theta_n1-phi_n1));

phi_n   = 2*pi.*rand([1,1]);

q_n = (2*C3*cos(phi_n)*sin(phi_n) + 2*C4*cos(phi_n)^2 - C4)/(2*C2*cos(phi_n) + 2*C1*sin(phi_n));
k_n = (2*C1*cos(phi_n)*q_n + 2*C2*sin(phi_n)*q_n - C3)/(4*q_n^2*cos(phi_n)^2-2*q_n^2);

x = -C1+cos(phi_n)*q_n*k_n;
y = -C2-sin(phi_n)*q_n*k_n;

r_n     = sqrt(x^2+y^2)/k_n;
theta_n = atan2(y,x);




check1 = C1 + k_n*r_n*cos(theta_n) - k_n*q_n*cos(phi_n)
check2 = C2 + k_n*r_n*sin(theta_n) + k_n*q_n*sin(phi_n)
check3 = C3 + 2*k_n*r_n*q_n*cos(theta_n-phi_n)
check4 = C4 + 2*k_n*r_n*q_n*sin(theta_n-phi_n)













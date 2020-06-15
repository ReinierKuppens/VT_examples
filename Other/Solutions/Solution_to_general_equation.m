
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

theta_n1 = 2*pi.*rand([1,nSpring-1]);
phi_n1   = 2*pi.*rand([1,nSpring-1]);

C1 = sum(r_n1.*cos(theta_n1))
C2 = sum(r_n1.*sin(theta_n1))

C3 = sum(q_n1.*cos(phi_n1))
C4 = sum(q_n1.*sin(phi_n1))

C5 = sum(q_n1.*r_n1.*cos(theta_n1-phi_n1))
C6 = sum(q_n1.*r_n1.*sin(theta_n1-phi_n1))




syms C1 C11 C2 C22 C3 C33 C4 C44 C5 C55 C6 C66
syms kn1 kn rn qn thetan phin 

Eq1 = C1 + C11*kn1 + kn*rn*cos(thetan);
Eq2 = C2 + C22*kn1 + kn*rn*sin(thetan);
Eq3 = C3 + C33*kn1 + kn*qn*cos(phin);
Eq4 = C4 + C44*kn1 + kn*qn*sin(phin);
Eq5 = C5 + C55*kn1 + kn*qn*rn*cos(thetan - phin);
Eq6 = C6 + C66*kn1 + kn*qn*rn*sin(thetan - phin);

syst = [Eq1;Eq2;Eq3;Eq4;Eq5;Eq6]

J = jacobian(syst,[kn1,kn,rn,qn,thetan,phin])





E1 = C1     +C11     +kn*rn*cos(thetan)
E2 = C2     +C22     +kn*rn*sin(thetan)
E3 = C3     +C33*qn1 +kn*qn*cos(phin)
E4 = C4     +C44*qn1 +kn*qn*sin(phin)
E5 = C5     +C55*qn1 +kn*rn*qn*cos(thetan-phin)
E6 = C6     +C66*qn1 +kn*rn*qn*sin(thetan-phin)






























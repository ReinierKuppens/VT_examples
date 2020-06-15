clear all 
close all 
clc 


nSpring  = 5;

rout     = 2;
rin      = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Sove the equation for coupling: sum krq*[cos(theta-phi);sin(theta-phi)]

%%%%%%%%%%%%%% 1

rbup = rin;

q_n2 =  2*rbup.*rand([1,nSpring-2]);
r_n2 =  2*rbup.*rand([1,nSpring-2]);
k_n2 =  2*rbup.*rand([1,nSpring-2]);

theta_n2 = 2*pi.*rand([1,nSpring-2]);
phi_n2   = 2*pi.*rand([1,nSpring-2]);



C1 = sum(r_n2.*k_n2.*cos(theta_n2));
C2 = sum(r_n2.*k_n2.*sin(theta_n2));

C3 = sum(q_n2.*k_n2.*cos(phi_n2));
C4 = sum(q_n2.*k_n2.*sin(phi_n2));

C5 = sum(q_n2.*r_n2.*k_n2.*cos(theta_n2-phi_n2));
C6 = sum(q_n2.*r_n2.*k_n2.*sin(theta_n2-phi_n2));

r_n1        =  2*rbup.*rand([1,1]);
k_n1        =  2*rbup.*rand([1,1]);
theta_n1    =  2*pi.*rand([1,1]);
phi_n1      =  2*pi.*rand([1,1]);

% C1 = 0
% C2 = 0
% C3 = 0
% C4 = 0
% C5 = 0
% C6 = 0


C11 = r_n1*k_n1*cos(theta_n1);
C22 = r_n1*k_n1*sin(theta_n1);
C33 = k_n1*cos(phi_n1);
C44 = k_n1*sin(phi_n1);
C55 = r_n1*k_n1*cos(theta_n1-phi_n1);
C66 = r_n1*k_n1*sin(theta_n1-phi_n1);


k_n = (-(((-C2-C22)*C5+C6*(C1+C11))^2*C33^2+(2*((C1+C11)*C5+C6*(C2+C22))*((-C2-C22)*C5+C6*(C1+C11))*C44+(2*(C2+C22)*((-C2-C22)*C5+C6*(C1+C11))*C3+2*C4*((C2+C22)*(C1+C11)*C5+C6*(C1^2+2*C1*C11+C11^2+2*C2^2+4*C2*C22+2*C22^2)))*C55-2*C66*((C1+C11)*((-C2-C22)*C5+C6*(C1+C11))*C3+2*C4*((C1^2+2*C1*C11+1/2*C2^2+C2*C22+C11^2+1/2*C22^2)*C5+1/2*C6*(C2+C22)*(C1+C11))))*C33+((C1+C11)*C5+C6*(C2+C22))^2*C44^2+(((2*(C2+C22)*(C1+C11)*C5-4*C6*(C1^2+2*C1*C11+1/2*C2^2+C2*C22+C11^2+1/2*C22^2))*C3-2*(C1+C11)*((C1+C11)*C5+C6*(C2+C22))*C4)*C55+2*C66*(((C1^2+2*C1*C11+C11^2+2*C2^2+4*C2*C22+2*C22^2)*C5-C6*(C2+C22)*(C1+C11))*C3-((C1+C11)*C5+C6*(C2+C22))*(C2+C22)*C4))*C44+(((-C2-C22)*C3+C4*(C1+C11))*C55+C66*((C1+C11)*C3+C4*(C2+C22)))^2)^(1/2)+(C1*C4+C11*C4-C2*C3-C22*C3)*C55+C66*(C1*C3+C11*C3+C2*C4+C22*C4)+(-C1*C44-C11*C44+C2*C33+C22*C33)*C5+(-C1*C33-C11*C33-C2*C44-C22*C44)*C6)/(-2*C5*C66+2*C55*C6) 



q_n1 = (-C1*C4-C11*C4+C2*C3+C22*C3+C6*k_n)/(C1*C44 + C11*C44 -C2*C33 - C22*C33 - C66*k_n)


x1 = -(C1+C11) / k_n
y1 = -(C2+C22) / k_n

x2 = -(C3 + C33*q_n1) / k_n
y2 = -(C4 + C44*q_n1) / k_n


theta_n = atan2(y1,x1);
phi_n   = atan2(y2,x2);

r_n = sqrt(x1^2+y1^2);
q_n = sqrt(x2^2+y2^2);

q = [q_n2,q_n1,q_n]
r = [r_n2,r_n1,r_n]
k = [k_n2,k_n1,k_n]

theta = [theta_n2,theta_n1,theta_n]
phi   = [phi_n2,phi_n1,phi_n]


check1 = sum(r.*k.*cos(theta))
check2 = sum(r.*k.*sin(theta))
check3 = sum(q.*k.*cos(phi))
check4 = sum(q.*k.*sin(phi))
check5 = sum(q.*r.*k.*cos(theta-phi))
check6 = sum(q.*r.*k.*sin(theta-phi))




















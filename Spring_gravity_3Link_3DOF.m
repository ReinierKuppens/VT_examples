


clear all 
close all 
clc 

alpha   = linspace(0,2*pi,150);

rho_1       = 1; 
rho_2       = -2;
rho_3       = 2;


rho_1       = 1; 
rho_2       = -5;
rho_3       = 2;

t1 = 0;
t2 = 1;
t3 = 1; 


nSpring03 = 2;
nSpring12 = 3 

k03     = 0.8.*rand([1,nSpring03-1]);
r03     = rand([1,nSpring03-1]);
q03     = rand([1,nSpring03-1]);
theta03 = 2*pi*rand([1,nSpring03-1]);
phi03   = 2*pi*rand([1,nSpring03-1]);

k03n     = 0.8.*rand;


n2 = 2*pi*rand;
v2 = rand;
m1 = 0.1.*rand;
m2 = 0.1.*rand;
m3 = 0.1.*rand;
g  = 9.81;



k12     = rand([1,nSpring12-1]); 
q12     = rand([1,nSpring12-1]); 
r12     = rand([1,nSpring12-1]); 
theta12 = 2*pi*rand([1,nSpring12-1]); 
phi12   = 2*pi*rand([1,nSpring12-1]); 

k12n        = 5.*rand;
r12n        = rand;
theta12n    = 2*pi*rand;


%% Theta03n and r03n 

C1 = sum(k03.*r03.*cos(theta03));
C2 = sum(k03.*r03.*sin(theta03));

xn = -C1;
yn = -C2;

theta03n = atan2(yn,xn);
r03n     = sqrt(xn^2+yn^2)/k03n;

check1 = C1+k03n*r03n*cos(theta03n);
check2 = C2+k03n*r03n*sin(theta03n);


%% phi03n and q03n

C3 = sum(k03.*q03.*t3.*cos(phi03))-g*m2*v2*cos(n2)-g*m3*t3;
C4 = sum(k03.*q03.*t3.*sin(phi03))+g*m2*v2*sin(n2) ;

xn = -C3;
yn = -C4;

phi03n = atan2(yn,xn);
q03n   = sqrt(xn^2+yn^2)/(k03n.*t3);

check3 = C3 + k03n.*q03n.*t3.*cos(phi03n);
check4 = C4 + k03n.*q03n.*t3.*sin(phi03n);

 

%%  n1 and v1 

C5 = sum(k03.*q03.*t2.*cos(phi03))+ k03n*q03n*t2*cos(phi03n)-g*m2*t2-g*m3*t2;
C6 = sum(k03.*q03.*t2.*sin(phi03))+ k03n*q03n*t2*sin(phi03n);

xn = C5;
yn = -C6;

n1 = atan2(yn,xn);
v1 = sqrt(xn^2+yn^2)/(g*m1);

check5 = C5 - (g*m1*v1)*cos(n1);
check6 = C6 + (g*m1*v1)*sin(n1);


%% n3 and v3

C7 = sum(k03.*r03.*q03.*cos(theta03-phi03)) + k03n*r03n*q03n*cos(theta03n-phi03n);
C8 = sum(k03.*r03.*q03.*sin(theta03-phi03)) + k03n*r03n*q03n*sin(theta03n-phi03n);

xn = C7;
yn = C8;

n3 = atan2(C8,C7);
v3 = sqrt(C8^2+C7^2)/(g*m3);

check7 = C7 - g*m3*v3*cos(n3);
check8 = C8 - g*m3*v3*sin(n3);

%% phi12n and q12n



C9  = -t2*t3*sum([k03,k03n]) +  sum(k12.*q12.*r12.*cos(theta12-phi12)-k12.*r12.*t2.*cos(theta12)) - k12n.*r12n.*t2.*cos(theta12n);
C10 =                           sum(k12.*q12.*r12.*sin(theta12-phi12)-k12.*r12.*t2.*sin(theta12)) - k12n.*r12n.*t2.*sin(theta12n);

xn = -C9;
yn = -C10;

gamma12n = atan2(yn,xn);
q12n     = sqrt(xn^2+yn^2)/(r12n*k12n);

phi12n   = -gamma12n + theta12n;

check9  = C9 + k12n*q12n*r12n*cos(theta12n-phi12n);
check10 = C10 + k12n*q12n*r12n*sin(theta12n-phi12n);


%% Make all points

q03         = [q03 q03n];
r03         = [r03 r03n];
theta03     = [theta03 theta03n];
phi03       = [phi03 phi03n];
k03         = [k03 k03n];

q12         = [q12 q12n];
r12         = [r12 r12n];
theta12     = [theta12 theta12n];
phi12       = [phi12 phi12n];
k12         = [k12 k12n];

%% Chech all equations

%rho2 + rho3

E1 = sum(k03.*r03.*t2.*cos(theta03));
E2 = sum(k03.*r03.*t2.*sin(theta03));

%rho3
E3 = sum(k03.*r03.*t3.*cos(theta03));
E4 = sum(k03.*r03.*t3.*sin(theta03));

%rho1 + rho2
E5 = sum(k03.*q03.*t3.*cos(phi03))-g*m2*v2*cos(n2)-g*m3*t3;
E6 = sum(k03.*q03.*t3.*sin(phi03))+g*m2*v2*sin(n2);

%rho1
E7 = sum(k03.*q03.*t2.*cos(phi03)) - g*m1*v1*cos(n1) - g*m2*t2 - g*m3*t2;
E8 = sum(k03.*q03.*t2.*sin(phi03)) + g*m1*v1*sin(n1) ;

%rho1 + rho2 +rho3
E9  = sum(k03.*r03.*q03.*cos(theta03-phi03)) - g*m3*v3*cos(n3);
E10 = sum(k03.*r03.*q03.*sin(theta03-phi03)) - g*m3*v3*sin(n3) ;

%rho2 
E11 = -sum(k03.*t2.*t3) + ( sum(k12.*q12.*r12.*cos(theta12-phi12) - k12.*r12.*t2.*cos(theta12)));
E12 =                       sum(k12.*q12.*r12.*sin(theta12-phi12) - k12.*r12.*t2.*sin(theta12)); 

CHECK = [E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12]
% keyboard 

%%


P03 = r03.*[cos(theta03);sin(theta03)];        
S03 = q03.*[cos(phi03);sin(phi03)]; 

P12 = r12.*[cos(theta12);sin(theta12)];
S12 = q12.*[cos(phi12);sin(phi12)];

T = @(a,tx)[cos(a),-sin(a),tx;
            sin(a), cos(a),0 ;
            0     , 0     ,1 ];
        
P03_prime         = cell(1,numel(alpha));
S03_prime         = cell(1,numel(alpha));        
P12_prime         = cell(1,numel(alpha));
S12_prime         = cell(1,numel(alpha));        

for k = 1:numel(alpha)
    P03_prime{k}          = T(rho_1*alpha(k),t1)*T(rho_2*alpha(k),t2)*T(rho_3*alpha(k),t3)*[P03;ones(1,length(P03(1,:)))];
    S03_prime{k}          = [S03;ones(1,length(S03(1,:)))];
    
    P12_prime{k}          = T(rho_1*alpha(k),t1)*T(rho_2*alpha(k),t2)*[P12;ones(1,length(P12(1,:)))];
    S12_prime{k}          = T(rho_1*alpha(k),t1)*[S12;ones(1,length(S12(1,:)))];
end



Com1 = [v1*cos(n1); 
        v1*sin(n1)];
Com2 = [v2*cos(n2); 
        v2*sin(n2)];
Com3 = [v3*cos(n3); 
        v3*sin(n3)];
    
C1_p = T(rho_1*alpha(1),t1)*[Com1;1];
C2_p = T(rho_1*alpha(1),t1)*T(rho_2*alpha(1),t2)*[Com2;1];
C3_p = T(rho_1*alpha(1),t1)*T(rho_2*alpha(1),t2)*T(rho_3*alpha(1),t3)*[Com3;1];




DNA           = initializeDNA;
DNA.incstr    = [1 3 6 3.*ones(1,nSpring12) 4.*ones(1,nSpring03)];
DNA.edgelabel = [1 1 1 2.*ones(1,nSpring12) 2.*ones(1,nSpring03)];


DNA.Mpar(1,:) = NaN;
DNA.Mpar(2,:) = [C1_p(1:2).' m1];
DNA.Mpar(3,:) = [C2_p(1:2).' m2];
DNA.Mpar(4,:) = [C3_p(1:2).' m3];


DNA.Hpar = [[t1;0],[t2+t1;0],[t3+t2+t1;0]];


DNA.Spar = [S12_prime{1}(1:2,:), S03_prime{1}(1:2,:);
            P12_prime{1}(1:2,:), P03_prime{1}(1:2,:);
            zeros(1,nSpring12+nSpring03);
            k12 k03];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                   make state            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


C1_prime         = cell(1,numel(alpha));
C2_prime         = cell(1,numel(alpha));
C3_prime         = cell(1,numel(alpha));

%  

for k = 1:numel(alpha)
    C1_prime{k}         = T(rho_1*alpha(k),t1)*[Com1;1];
    C2_prime{k}         = T(rho_1*alpha(k),t1)*T(rho_2*alpha(k),t2)*[Com2;1];
    C3_prime{k}         = T(rho_1*alpha(k),t1)*T(rho_2*alpha(k),t2)*T(rho_3*alpha(k),t3)*[Com3;1];
end

% keyboard 

for k = 1:length(alpha)
    state(k,:) = [  C1_prime{k}(1:2,:).', (rho_1)*alpha(k), ...
                    C2_prime{k}(1:2,:).', (rho_1 + rho_2)*alpha(k), ...
                    C3_prime{k}(1:2,:).', (rho_1 + rho_2 + rho_3)*alpha(k)];
end
t=alpha;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                   test    energies             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for k = 1:numel(alpha) 
    for ii = 1:nSpring03
        U03i(k,ii) = 0.5*k03(ii)*((P03_prime{k}(1,ii)-S03_prime{k}(1,ii))^2 + (P03_prime{k}(2,ii)-S03_prime{k}(2,ii))^2);
    end
    for ii = 1:nSpring12
        U12i(k,ii) = 0.5*k12(ii)*((P12_prime{k}(1,ii)-S12_prime{k}(1,ii))^2 + (P12_prime{k}(2,ii)-S12_prime{k}(2,ii))^2);
    end
end

Utotal = sum(U03i,2) + sum(U12i,2);

Ug     = getGravityEnergy(DNA,t,state).';

% keyboard 
Utotal = Utotal + sum(Ug,2); 

% 
% figure;hold on 
% for k = 1:nSpring03
%     plot(alpha,U03i(:,k),'r','linewidth',1)
% end

% for k = 1:nSpring12
%     plot(alpha,U12i(:,k),'y','linewidth',1)
% end
% plot(alpha,Ug,'--k','linewidth',1)
% plot(alpha,Utotal,'-k','linewidth',1)
% title('test1')

% keyboard 


[Es,dEs]            = getEnergies(DNA,t,state);
% figure;hold on;
Ug = Ug.';
Et = zeros(1,numel(t));
for k = 1:numel(Es)
%     plot(alpha,Es{k},'r')
    Et = Et + Es{k};
end

Et = Et+sum(Ug);
% plot(alpha,Et,'k')
% plot(alpha,Ug,'--k')
% title('test 2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                     ANIMATE                    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

warning off 

plotEnergies(DNA,t,state)
plotmDNA(DNA,t,state)

% animateDNA(DNA,t,state)
% 

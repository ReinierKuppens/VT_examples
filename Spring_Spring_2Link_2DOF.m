
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           Parameter Definitions                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSpring = 5;

r1      = 2;
r2      = 1;

rho     = -r1/r2;
% rho     = 11;


alpha   = linspace(0,2*pi,150);
beta    = rho*alpha;


%============== One solution with nS=3 ====================================
nSpring = 3;

r3 = 1;
r4 = 3;
stiffness = ones(1,nSpring) ;


% theta = linspace(0,2*pi,nSpring+1)
% theta = theta(1:nSpring) 
% theta = [theta(1) theta(end) theta(2:end-1)]

theta = [(3/12)*pi (19/12)*pi (11/12)*pi];

% phi = linspace(0,2*pi,nSpring+1) 
% phi = phi(1:nSpring)
% phi = [phi(1) phi(1:end-1)]

phi   = [(3/24)*pi  (19/24)*pi (35/24)*pi]; 

P1 = r3.*[cos(theta) ;sin(theta)];        
P2 = r4.*[cos(phi)   ;sin(phi)]; 
%==========================================================================


 

circle1 = [r1*cos(0:0.01:2*pi)        ;r1*sin(0:0.01:2*pi)];
circle2 = [r2*cos(0:0.01:2*pi)        ;r2*sin(0:0.01:2*pi)];
circle3 = [(r1-r2)*cos(0:0.01:2*pi)   ;(r1-r2)*sin(0:0.01:2*pi)];

T = @(a,tx)[cos(a),-sin(a),tx;
            sin(a), cos(a),0 ;
            0     , 0     ,0 ];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%     Compute  distance^2 and energy             %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
circle2_prime   = cell(1,numel(alpha));
P_prime         = cell(1,numel(alpha));
%Compute point locations of complete domain alpha
for k = 1:numel(alpha)
    
    circle2_prime{k}    = T(alpha(k),0)*T(beta(k),(r1-r2))*[circle2;ones(1,length(circle2(1,:)))];
    P1_prime{k}          = T(alpha(k),0)*T(beta(k),(r1-r2))*[P1;ones(1,length(P1(1,:)))];
    P2_prime{k}          = [P2;ones(1,length(P1(1,:)))];

    %     circle2_prime{k}    = R(-alpha(k)*rho)*circle2 + [(r1-r2)*cos(alpha(k));(r1-r2)*sin(alpha(k))];
%     P_prime{k}          = R(-alpha(k)*rho)*P       + [(r1-r2)*cos(alpha(k));(r1-r2)*sin(alpha(k))];
end

%Compute energies
% stifness = ones(1,nSpring) ;
for k = 1:numel(alpha)
    for ii = 1:nSpring
        Ui(k,ii) = 0.5*stiffness(ii)*((P1_prime{k}(1,ii)-P2_prime{k}(1,ii))^2 + (P1_prime{k}(2,ii)-P2_prime{k}(2,ii))^2);
    end
end
Utotal = sum(Ui,2);


% figure;hold on; 
% plot(alpha,Ui(:,:)); plot(alpha,Utotal,'--r')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                Plot Mechanism                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                     Make DNA                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DNA           = initializeDNA;
DNA.incstr    = [1 3 2.*ones(1,nSpring)];
DNA.edgelabel = [1 1 2.*ones(1,nSpring)];

DNA.Mpar(1,:) = NaN;
DNA.Mpar(2,:) = [(r1-r2)/2 0 0];
DNA.Mpar(3,:) = [(r1-r2)   0 0];

DNA.Hpar = [[0;0],[r1-r2;0]];

DNA.Spar = [P2_prime{1}(1:2,:);
    P1_prime{1}(1:2,:);
    ones(1,nSpring);
    stiffness];

[~,comPoint] = getMassLocations(DNA);

m1 = comPoint{1}.';
m2 = comPoint{2}.'-[r1-r2;0];

M1_prime         = cell(1,numel(alpha));
M2_prime         = cell(1,numel(alpha));
for k = 1:numel(alpha)
    M1_prime{k}         = T(alpha(k),0)*[m1;1];
    M2_prime{k}         = T(alpha(k),0)*T(alpha(k)*rho,(r1-r2))*[m2;1];
end

for k = 1:length(alpha)
    state(k,:) = [M1_prime{k}(1:2,:).', alpha(k), M2_prime{k}(1:2,:).', (rho+1)*alpha(k) ];
%     t(k) = k;
end
t=alpha;

warning off
figure('position',[2200,200,350,350])
plotmDNA(DNA,t,state)
% keyboard
% 

[Es,dEs]            = getEnergies(DNA,t,state);
% keyboard
figure;hold on;
% keyboard

Et = zeros(1,numel(t));
for k = 1:numel(Es)
    plot(Es{k},'r')
    Et = Et + Es{k};
end
plot(Et,'k')


plotTrajectory(DNA,t,state)

plotEnergies(DNA,t,state) 

animateDNA(DNA,t,state)









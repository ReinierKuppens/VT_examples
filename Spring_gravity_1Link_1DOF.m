
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           Parameter Definitions                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSpring = 3;
alpha   = linspace(0,2*pi,150); 
g       = 9.81;

lb = 0; 
ub = 2; 

axi = lb + (ub-lb).*rand(1,nSpring);
ayi = lb + (ub-lb).*rand(1,nSpring);
bx  =  sum(axi);
by  =  sum(ayi);

ai      = sqrt(axi.^2 + ayi.^2);
gammai  = atan2(ayi,axi);

ki = 1 + (2-1).*rand(1,nSpring);
ri = 0 + (2-0).*rand(1,nSpring);
qi = ai./(ki.*ri);

phii    = 2*pi.*rand(1,nSpring);
thetai  = gammai + phii;  

b       = sqrt(bx.^2 + by.^2);
eta1    = atan2(by,bx);

m1      = 0.1 + (0.5-0.1).*rand(1); 
v1      = b/(g*m1);

P = ri.*[cos(thetai);sin(thetai)];        
S = qi.*[cos(phii);sin(phii)]; 

T = @(a,tx)[cos(a),-sin(a),tx;
            sin(a), cos(a),0 ;
            0     , 0     ,0 ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%     Compute  distance^2 and energy             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
P_prime         = cell(1,numel(alpha));
S_prime         = cell(1,numel(alpha));

%Compute point locations of complete domain alpha
for k = 1:numel(alpha)
    P_prime{k}          = T(alpha(k),0)*[P;ones(1,length(P(1,:)))];
    S_prime{k}          = [S;ones(1,length(S(1,:)))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%     Construct Mechanism DNA datastructure      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DNA           = initializeDNA;
DNA.incstr    = [1  1.*ones(1,nSpring)];
DNA.edgelabel = [1  2.*ones(1,nSpring)];
DNA.Mpar(1,:) = NaN;
DNA.Mpar(2,:) = [v1*cos(eta1), v1*sin(eta1) m1];
DNA.Hpar      = [0;0];
DNA.Spar       = [S_prime{1}(1:2,:);
                  P_prime{1}(1:2,:);
                  zeros(1,nSpring);
                  ki];
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                   make state            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        
[~,comPoint] = getMassLocations(DNA);
C1           = comPoint{1}.';
C1_prime     = cell(1,numel(alpha));

for k = 1:numel(alpha)
    C1_prime{k}         = T(alpha(k),0)*[C1;1];
end

for k = 1:length(alpha)
    state(k,:) = [C1_prime{k}(1:2,:).', alpha(k)];
end
t=alpha;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                     ANIMATE                    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off 
plotEnergies(DNA,t,state)
plotmDNA(DNA,t,state)
animateDNA(DNA,t,state)









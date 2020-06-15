
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           Parameter Definitions                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = linspace(0,2*pi,150); 

nSpring = 3;

g = 9.81;

lb = 0; 
ub = 2; 

axi = lb + (ub-lb).*rand(1,nSpring);
ayi = lb + (ub-lb).*rand(1,nSpring);

bxi = lb + (ub-lb).*rand(1,nSpring);
byi = lb + (ub-lb).*rand(1,nSpring);

cx =  lb + (ub-lb).*rand(1,1);
cy = sum(ayi) + sum(byi);

d = sum(axi) - sum(bxi) - cx;


sum(axi)-sum(bxi)-cx -d;
sum(ayi)+sum(byi)-cy;

% keyboard 


ai = sqrt(axi.^2 + ayi.^2);
bi = sqrt(bxi.^2 + byi.^2);

gammai = atan2(ayi,axi);
thetai = atan2(byi,bxi);

c = sqrt(cx.^2+cy.^2);
eta1 = atan2(cy,cx);


phii = gammai + thetai; 

m2 = 0.1 + (0.5-0.1).*rand(1,1);
t2 = d/(m2*g); 
 
m1 = lb + (ub-lb).*rand(1,1);
v1 = c/(g*m1);

ki = 1 + (10-1).*rand(1,nSpring);
ri = bi./(ki.*t2);
qi = ai./(ri.*ki);

eta2 = 2*pi.*rand(1,1);
v2   = lb + (ub-lb).*rand(1,1);



P = ri.*[cos(thetai);sin(thetai)];        
S = qi.*[cos(phii);sin(phii)]; 

% circle1 = [r1*cos(0:0.01:2*pi)        ;r1*sin(0:0.01:2*pi)];
% circle2 = [r2*cos(0:0.01:2*pi)        ;r2*sin(0:0.01:2*pi)];
% circle3 = [(r1-r2)*cos(0:0.01:2*pi)   ;(r1-r2)*sin(0:0.01:2*pi)];

T = @(a,tx)[cos(a),-sin(a),tx;
            sin(a), cos(a),0 ;
            0     , 0     ,0 ];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%     Compute  distance^2 and energy             %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% circle2_prime   = cell(1,numel(alpha));
P_prime         = cell(1,numel(alpha));
S_prime         = cell(1,numel(alpha));

rho_1 = 1; 
t1    = 0;
rho_2 = -1;

%Compute point locations of complete domain alpha
for k = 1:numel(alpha)
    P_prime{k}          = T(rho_1*alpha(k),t1)*T(rho_2*alpha(k),t2)*[P;ones(1,length(P(1,:)))];
    S_prime{k}          = T(rho_1*alpha(k),t1)*[S;ones(1,length(S(1,:)))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                Plot Mechanism                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                     Make DNA                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DNA           = initializeDNA;
DNA.incstr    = [1 3 3.*ones(1,nSpring)];
DNA.edgelabel = [1 1 2.*ones(1,nSpring)];


C1 = [v1*cos(eta1); v1*sin(eta1)]
C2 = [v2*cos(eta2); v2*sin(eta2)]

C1_p = T(rho_1*alpha(1),t1)*[C1;1];
C2_p = T(rho_1*alpha(1),t1)*T(rho_2*alpha(1),t2)*[C2;1]

% keyboard 

DNA.Mpar(1,:) = NaN;
DNA.Mpar(2,:) = [C1_p(1:2).' m1];
DNA.Mpar(3,:) = [C2_p(1:2).' m2];




DNA.Hpar = [[t1;0],[t2;0]];

DNA.Spar = [S_prime{1}(1:2,:);
            P_prime{1}(1:2,:);
            zeros(1,nSpring);
            ki];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                   make state            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        
% [~,comPoint] = getMassLocations(DNA);
% 
% C1 = comPoint{1}.';
% C2 = comPoint{2}.';
% keyboard 

C1_prime         = cell(1,numel(alpha));
C2_prime         = cell(1,numel(alpha));

for k = 1:numel(alpha)
    C1_prime{k}         = T(rho_1*alpha(k),t1)*[C1;1];
    C2_prime{k}         = T(rho_1*alpha(k),t1)*T(rho_2*alpha(k),t2)*[C2;1];
end

for k = 1:length(alpha)
    state(k,:) = [C1_prime{k}(1:2,:).', rho_1*alpha(k), C2_prime{k}(1:2,:).', (rho_2+rho_1)*alpha(k)];
%     t(k) = k;
end
t=alpha;

% warning off
% figure('position',[2200,200,350,350])
% plotmDNA(DNA,t,state)
% hold on

% plot(state(:,1),state(:,2))
% hold on; plot(state(:,4),state(:,5))
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                   test    energies             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %Compute energies
% stifness = ones(1,nSpring) ;
for k = 1:numel(alpha)
    for ii = 1:nSpring
        Ui(k,ii) = 0.5*ki(ii)*((P_prime{k}(1,ii)-S_prime{k}(1,ii))^2 + (P_prime{k}(2,ii)-S_prime{k}(2,ii))^2);
    end
end
Utotal = sum(Ui,2);

Ug     = getGravityEnergy(DNA,t,state).';
Utotal = Utotal + sum(Ug,2); 

% figure;hold on 
% for k = 1:nSpring
%    plot(alpha,Ui(:,k),'r','linewidth',1) 
% end
% plot(alpha,Ug,'--k','linewidth',1)
% plot(alpha,Utotal,'-k','linewidth',1)
% title('test1')

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

figure('color',[1,1,1]); 
plotmDNA(DNA,t,state)

animateDNA(DNA,t,state)









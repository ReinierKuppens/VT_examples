


clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           Parameter Definitions                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSpring = 2;

alpha   = linspace(0,2*pi,150);


rho_2       = 2;
rho_1       = 1; 
t1          = 0;

%Link parameters, randomly defined. 
g       = 9.81;

m1      = 0.1 +(0.5-0.1)*rand([1,1]);
eta1    = 2*pi*rand([1,1]);
v1      = 0.1 +(2-0.1)*rand([1,1]);

m2   = 0.1 +(0.5-0.1)*rand([1,1]);
eta2 = 2*pi*rand([1,1]);
v2 = 0.1 +(2-0.1)*rand([1,1]);

rbup = 1;

q_n1 =  0.1+2*rbup.*rand([1,nSpring-1]);
r_n1 =  0.1+2*rbup.*rand([1,nSpring-1]);
k_n1 =  rbup +  2*rbup.*rand([1,nSpring-1]);

theta_n1 = 2*pi.*rand([1,nSpring-1]);
phi_n1   = 2*pi.*rand([1,nSpring-1]);


% Define Coefficients
C1 = sum(r_n1.*k_n1.*cos(theta_n1));
C2 = sum(r_n1.*k_n1.*sin(theta_n1));

C3 = sum(q_n1.*k_n1.*cos(phi_n1))-g*m2;
C4 = sum(q_n1.*k_n1.*sin(phi_n1));

C5 = sum(q_n1.*r_n1.*k_n1.*cos(theta_n1-phi_n1)) - g*m2*v2*cos(eta2);
C6 = sum(q_n1.*r_n1.*k_n1.*sin(theta_n1-phi_n1)) - g*m2*v2*sin(eta2);


%% We can do: 
k_n = -(C1^2+C2^2)*(C3*sin(eta1)+C4*cos(eta1))/(C1*C5*sin(eta1)-C1*C6*cos(eta1)+C2*C5*cos(eta1)+C2*C6*sin(eta1));
t2  = -g*m1*v1*(C1*sin(eta1)+C2*cos(eta1))/(C1*C4-C2*C3-C6*k_n);



%% Or we do: This is preferred because k>=0; 
k_n  = 5; 
eta1 = -atan2((C1^2*C4-C1*C6*k_n+C2^2*C4+C2*C5*k_n),(C1^2*C3+C1*C5*k_n+C2^2*C3+C2*C6*k_n));
t2   = -g*m1*v1*(C1*sin(eta1)+C2*cos(eta1))/(C1*C4-C2*C3-C6*k_n);


%%

x1 =  -C1/k_n;
y1 =  -C2/k_n;

x2 =   (-C3+g*m1*v1/t2*cos(eta1))/k_n;
y2 =   (-C4-g*m1*v1/t2*sin(eta1))/k_n;

theta_n = atan2(y1,x1);
phi_n   = atan2(y2,x2);
r_n     = sqrt(x1^2+y1^2);
q_n     = sqrt(x2^2+y2^2);

% theta_n = atan2(C2,C1);
% r_n     = -C1/(k_n*cos(theta_n))

% phi_n   = 

 

qi = [q_n1,q_n];
ri = [r_n1,r_n];
ki = [k_n1,k_n];

thetai = [theta_n1,theta_n];
phii   = [phi_n1,phi_n];


check1 = sum(ri.*ki.*cos(thetai)) ;
check2 = sum(ri.*ki.*sin(thetai)) ;
check3 = sum(qi.*ki.*cos(phii)) - g*m2 - (1/t2)*g*m1*v1*cos(eta1);
check4 = sum(qi.*ki.*sin(phii)) + (1/t2)*g*m1*v1*sin(eta1);
check5 = sum(qi.*ri.*ki.*cos(thetai-phii)) - g*m2*v2*cos(eta2);
check6 = sum(qi.*ri.*ki.*sin(thetai-phii)) - g*m2*v2*sin(eta2);

check = [check1, check2, check3, check4, check5, check6]


% keyboard 


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



%Compute point locations of complete domain alpha
for k = 1:numel(alpha)
    P_prime{k}          = T(rho_1*alpha(k),t1)*T(rho_2*alpha(k),t2)*[P;ones(1,length(P(1,:)))];
    S_prime{k}          = [S;ones(1,length(S(1,:)))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                Plot Mechanism                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                     Make DNA                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DNA           = initializeDNA;
DNA.incstr    = [1 3 2.*ones(1,nSpring)];
DNA.edgelabel = [1 1 2.*ones(1,nSpring)];


C1 = [v1*cos(eta1); 
      v1*sin(eta1)];
C2 = [v2*cos(eta2); 
      v2*sin(eta2)];

C1_p = T(rho_1*alpha(1),t1)*[C1;1];
C2_p = T(rho_1*alpha(1),t1)*T(rho_2*alpha(1),t2)*[C2;1];

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
    state(k,:) = [  C1_prime{k}(1:2,:).', rho_1*alpha(k),... 
                    C2_prime{k}(1:2,:).', (rho_2+rho_1)*alpha(k)];
%     t(k) = k;
end
t=alpha;

% warning off
% figure('position',[2200,200,350,350])
% plotmDNA(DNA,t,state)
% hold on
% 
% plot(state(:,1),state(:,2))
% hold on; plot(state(:,4),state(:,5))
 
% keyboard 

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




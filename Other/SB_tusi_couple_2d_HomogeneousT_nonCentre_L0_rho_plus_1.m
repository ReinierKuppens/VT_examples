
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           Parameter Definitions                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSpring = 2;

r1      = 4;
r2      = 2;

rho     = 1;

alpha   = linspace(0,2*pi,200);
beta    = rho*linspace(0,2*pi,200);

rbup = r2;

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


q           = [q_n1 q_n];
r           = [r_n1 r_n];
stiffness   = [k_n1 k_n];

theta       = [theta_n1 theta_n];
phi         = [phi_n1 phi_n];


%==========================================================================


P1 = r.*[cos(theta);sin(theta)];        
P2 = q.*[cos(phi);sin(phi)];  

circle1 = [r1*cos(0:0.01:2*pi)        ;r1*sin(0:0.01:2*pi)];
circle2 = [r2*cos(0:0.01:2*pi)        ;r2*sin(0:0.01:2*pi)];
circle3 = [(r1-r2)*cos(0:0.01:2*pi)   ;(r1-r2)*sin(0:0.01:2*pi)];

T = @(a,tx)[cos(a),-sin(a),tx;
            sin(a), cos(a),0 ;
            0     , 0     ,0 ];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%     Compute  distance^2 and energy             %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
DNA.Mpar(2,:) = [(r1-r2)/2 0 1];
DNA.Mpar(3,:) = [0 0 1];

DNA.Hpar = [[0;0],[r1-r2;0]];

DNA.Spar = [P2_prime{1}(1:2,:);
    P1_prime{1}(1:2,:);
    zeros(1,nSpring);
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

% keyboard
animateDNA(DNA,t,state)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                Initial Plots                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
r3 = r;
r4 = q;

fh      = figure('color',[1,1,1],'position',[100,100,1000,500]);subplot(1,2,1);hold on ;axis equal;
ht      = title(['R1 = ',num2str(r1),', R2 = ',num2str(r2),', Offset = ',num2str(r3)]); %Print time in title

f       = getframe(gcf);
sizef   = size(f.cdata);
width   = sizef(1);
height  = sizef(2);
% range = max([r1,r2,2*r3]);
% axis([-range*1.2,range*1.2,-range*1.2,range*1.2])

subplot(1,2,1); hold on
H1 = plot(circle1(1,:),circle1(2,:),'k','Linewidth',1);
H3 = plot(circle2(1,:)+circle3(1,1),circle2(2,:)+circle3(2,1),'g');
% H4 = plot([P1(1,:),P2(1,:)]+circle3(1,1),[P1(2,:),P2(2,:)]+circle3(2,1),'b.','markersize',20);

H4 = plot(P1_prime{1}(1,:),P1_prime{1}(2,:),'g.','markersize',20);
H5 = plot(P2_prime{1}(1,:),P2_prime{1}(2,:),'b.','markersize',20);
H6 = plot([P1_prime{1}(1,:);P2_prime{1}(1,:)],[P1_prime{1}(2,:);P2_prime{1}(2,:)],'k');

% keyboard 

plot([0,0],[-r1,r1],'color',[0.7,0.7,0.7])
plot(circle3(1,:),circle3(2,:),'color',[0.7,0.7,0.7])
xlabel('Distance [m]')
ylabel('Distance [m]')
% keyboard 
for k = 1:numel(alpha)
    plot(P1_prime{k}(1,:),P1_prime{k}(2,:),'.r','markersize',5)
    plot(P2_prime{k}(1,:),P2_prime{k}(2,:),'.b','markersize',2)
end

subplot(1,2,2)

plot(abs(alpha)+abs(beta),Utotal);hold on;
for k = 1:length(P1_prime{1}(1,:))
    plot(abs(alpha)+abs(beta),Ui(:,k));
end

U1 = plot(alpha(1),Utotal(1),'.b','markersize',20);
for ii = 1:length(P1_prime{k}(1,:))
    U2(ii) = plot(abs(alpha(1))+abs(beta(1)),Ui(1,ii),'.r','markersize',20);
end
xlabel('angle [rad]')
ylabel('Energy [J]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                Make animation                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mCounter=1;
for k = 1:numel(alpha)
    
    set(H3,'Xdata',circle2_prime{k}(1,:),'Ydata',circle2_prime{k}(2,:))
    set(H4,'Xdata',P1_prime{k}(1,:),'Ydata',P1_prime{k}(2,:))
    set(H5,'Xdata',P2_prime{k}(1,:),'Ydata',P2_prime{k}(2,:))
    for jj = 1:numel(H6)
        set(H6(jj),'Xdata',[P1_prime{k}(1,jj);P2_prime{k}(1,jj)],'Ydata',[P1_prime{k}(2,jj);P2_prime{k}(2,jj)])
    end
    set(U1,'Xdata',abs(alpha(k))+abs(beta(k)),'Ydata',Utotal(k))
    
for ii = 1:length(P1_prime{k}(1,:))
    set(U2(ii),'Xdata',abs(alpha(k))+abs(beta(k)),'Ydata',Ui(k,ii))
end
    
    f=getframe(gcf);
    if mCounter == 1
        [mov(:,:,1,mCounter), map] = rgb2ind(f.cdata, 256, 'nodither');
    else
               mov(:,:,1,mCounter) = rgb2ind(f.cdata, map, 'nodither');
    end
    
    mCounter = mCounter+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Save animations                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

directory = 'C:\Users\Reinier\Box Sync\Reinier PhD\3_Matlab\SB_Tusi_couple\Animations';
directory = 'D:\pkuppens\Box Sync\Reinier PhD\3_Matlab\SB_Tusi_couple\Animations';

idx=0;
cd(directory)
existingNames=dir('*.gif');
for k = 1:numel(existingNames)
    idx(k)=str2num(existingNames(k).name(end-7:end-4));
end
maxIdx = max(idx);

animation_name = strcat(directory,'\Animation',num2str(maxIdx+1.','%04d'),'.gif');
imwrite(mov, map,animation_name, 'DelayTime', 0, 'LoopCount', inf);
% 
% 
% 
% 
% 











clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           Parameter Definitions                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSpring = 5;

rout     = 2;
rin      = 1;

rho     =  rout/rin;

alpha   = linspace(0,2*pi,200);
beta    = -rho*linspace(0,2*pi,200);


%============== random solution ====================================

% ============ 1
rbup = rin;

x1 = -rbup + (2*rbup).*rand([1,nSpring-1]);
y1 = -rbup + (2*rbup).*rand([1,nSpring-1]);

x1(end+1) = -sum(x1);
y1(end+1) = -sum(y1);

A1      = sqrt(x1.^2 + y1.^2)
theta   = atan2(y1,x1);

sum(A1.*cos(theta))
sum(A1.*sin(theta))

r   = rbup.*rand([1,nSpring]);

stiffness = A1./((rin-rout).*r)

r.*(rin-rout).*stiffness

% ============ 2

x2 = -rbup + (2*rbup).*rand([1,nSpring-1]);
y2 = -rbup + (2*rbup).*rand([1,nSpring-1]);

x2(end+1) = -sum(x2);
y2(end+1) = -sum(y2);

A2      = sqrt(x2.^2 + y2.^2)

sum(A2.*cos(theta))
sum(A2.*sin(theta))

rho2 = -rin/rout

q       = A2./(stiffness.*rho2.*(rin-rout))
phi     = atan2(y2,x2);

sum(rho2.*q.*stiffness.*(rin-rout).*cos(phi))
sum(rho2.*q.*stiffness.*(rin-rout).*sin(phi))

sum((rho2+1).*stiffness.*r.*q.*cos(theta-phi))
sum((rho2+1).*stiffness.*r.*q.*sin(theta-phi))


% keyboard 

%============== One solution with nS=3 ====================================
r3 = 1;
r4 = 2.5;
stiffness = ones(1,nSpring) ;
% 
% theta = [(3/12)*pi (19/12)*pi (11/12)*pi];
% phi   = [0          (8/12)*pi (16/12)*pi]+pi/8;
%==========================================================================


P1 = r3.*[cos(theta);sin(theta)];        
P2 = r4.*[cos(phi);sin(phi)];  

circle1 = [rout*cos(0:0.01:2*pi)        ;rout*sin(0:0.01:2*pi)];
circle2 = [rin*cos(0:0.01:2*pi)        ;rin*sin(0:0.01:2*pi)];
circle3 = [(rout-rin)*cos(0:0.01:2*pi)   ;(rout-rin)*sin(0:0.01:2*pi)];

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
    
    circle2_prime{k}    = T(alpha(k),0)*T(-beta(k),(rout-rin))*[circle2;ones(1,length(circle2(1,:)))];
    P1_prime{k}          = T(alpha(k),0)*T(-beta(k),(rout-rin))*[P1;ones(1,length(P1(1,:)))];
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
DNA.Mpar(2,:) = [(rout-rin)/2 0 0];
DNA.Mpar(3,:) = [(rout-rin)   0 0];

DNA.Hpar = [[0;0],[rout-rin;0]];

DNA.Spar = [P2_prime{1}(1:2,:);
    P1_prime{1}(1:2,:);
    2.*ones(1,nSpring);
    stiffness];

[~,comPoint] = getMassLocations(DNA);

m1 = comPoint{1}.';
m2 = comPoint{2}.'-[rout-rin;0];

M1_prime         = cell(1,numel(alpha));
M2_prime         = cell(1,numel(alpha));
for k = 1:numel(alpha)
    M1_prime{k}         = T(alpha(k),0)*[m1;1];
    M2_prime{k}         = T(alpha(k),0)*T(-alpha(k)*rho,(rout-rin))*[m2;1];
end

for k = 1:length(alpha)
    state(k,:) = [M1_prime{k}(1:2,:).', alpha(k), M2_prime{k}(1:2,:).', -(rho-1)*alpha(k) ];
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


fh      = figure('color',[1,1,1],'position',[100,100,1000,500]);subplot(1,2,1);hold on ;axis equal;
ht      = title(['R1 = ',num2str(rout),', R2 = ',num2str(rin),', Offset = ',num2str(r3)]); %Print time in title

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

plot([0,0],[-rout,rout],'color',[0.7,0.7,0.7])
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










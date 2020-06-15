
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           Parameter Definitions                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSpring     = 7;
r1          = 2;
r2          = 2;
r3          = 2;
rho         = -(r1/r2);
theta       = linspace(0,2*pi,nSpring+1);theta(end)=[];
P           = r3.*[cos(theta);sin(theta)];
stiffness   = ones(1,nSpring);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%     Evely spaced spring insertions          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSpring     = 5;
r1          = 2;
r2          = 1;
rho         = (r1/r2);

A           = ones(1,nSpring);
r           = r3.*rand([1,nSpring]);
r(1)        = r2;

theta       = linspace(0,2*pi,nSpring+1);theta(end)=[];
P           = r.*[cos(theta);sin(theta)];
stiffness   = A./r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%     Random spring insertions        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% x = -r3 + (2*r3).*rand([1,nSpring-1]);
% y = -r3 + (2*r3).*rand([1,nSpring-1]);
% r = r3.*rand([1,nSpring]);
% % r(1) = r3;
% 
% x(end+1) = -sum(x);
% y(end+1) = -sum(y);
% 
% A           = sqrt(x.^2 + y.^2);
% theta       = atan2(y,x);
% stiffness   = A./r;
% P           = r.*[cos(theta);sin(theta)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%             some more definitions              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha   = linspace(0,2*pi,200);
beta    = alpha.*rho;

circle1 = [r1*cos(alpha)        ;r1*sin(alpha)];
circle2 = [r2*cos(alpha)        ;r2*sin(alpha)];
circle3 = [(r1-r2)*cos(alpha)   ;(r1-r2)*sin(alpha)];

T = @(a,tx)[cos(a),-sin(a),tx;
    sin(a), cos(a),0 ;
    0     , 0     ,0 ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%     Compute  distance^2 and energy             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

circle2_prime    = cell(1,numel(alpha));
P_prime          = cell(1,numel(alpha));
%Compute point locations of complete domain alpha
for k = 1:numel(alpha)
    circle2_prime{k}    = T(alpha(k),0)*T(-alpha(k)*rho,(r1-r2))*[circle2;ones(1,length(circle2(1,:)))];
    P_prime{k}          = T(alpha(k),0)*T(-alpha(k)*rho,(r1-r2))*[P;ones(1,length(P(1,:)))];
end

%Compute energies
for k = 1:numel(beta)
    for ii = 1:nSpring
        Ui(k,ii) = 0.5*stiffness(ii)*(P_prime{k}(1,ii)^2 + P_prime{k}(2,ii)^2);
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

DNA.Spar = [zeros(2,nSpring);
            P_prime{1}(1:2,:);
            zeros(1,nSpring);
            stiffness];

DNA.rho = rho; 

[~,comPoint] = getMassLocations(DNA);

m1 = comPoint{1}.';
m2 = comPoint{2}.'-[r1-r2;0];

M1_prime         = cell(1,numel(alpha));
M2_prime         = cell(1,numel(alpha));
for k = 1:numel(alpha)
    M1_prime{k}         = T(alpha(k),0)*[m1;1];
    M2_prime{k}         = T(alpha(k),0)*T(-alpha(k)*rho,(r1-r2))*[m2;1];
end

for k = 1:length(alpha)
    state(k,:) = [M1_prime{k}(1:2,:).', alpha(k), M2_prime{k}(1:2,:).', -(rho-1)*alpha(k) ];
    t(k) = k;
end
t = alpha;

warning off
figure('position',[2200,200,350,350])
plotmDNA(DNA,t,state)

animateDNA(DNA,t,state)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                Initial Plots                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% fh      = figure('color',[1,1,1],'position',[100,100,1000,500]);subplot(1,2,1);hold on ;axis equal;
% ht      = title(['R1 = ',num2str(r1),', R2 = ',num2str(r2),', Offset = ',num2str(r3)]); %Print time in title
% 
% f       = getframe(gcf);
% sizef   = size(f.cdata);
% width   = sizef(1);
% height  = sizef(2);
% 
% subplot(1,2,1); hold on
% H1 = plot(circle1(1,:),circle1(2,:),'k','Linewidth',1);
% H3 = plot(circle2(1,:)+circle3(1,1),circle2(2,:)+circle3(2,1),'g');
% H4 = plot(P(1,:)+circle3(1,1),P(2,:)+circle3(2,1),'b.','markersize',20);
% 
% plot([0,0],[-r1,r1],'color',[0.7,0.7,0.7])
% plot(circle3(1,:),circle3(2,:),'color',[0.7,0.7,0.7])
% xlabel('Distance [m]')
% ylabel('Distance [m]')
% 
% for k = 1:numel(beta)
%     plot(P_prime{k}(1,:),P_prime{k}(2,:),'.r','markersize',5)
% end
% 
% subplot(1,2,2)
% 
% plot(alpha,Utotal);hold on;
% for k = 1:length(P_prime{1}(1,:))
%     plot(alpha,Ui(:,k));
% end
% 
% U1 = plot(alpha(1),Utotal(1),'.b','markersize',20);
% for ii = 1:length(P_prime{k}(1,:))
%     U2(ii) = plot(alpha(1),Ui(1,ii),'.r','markersize',20);
% end
% xlabel('angle [rad]')
% ylabel('Energy [J]')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                Make animation                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % mCounter=1;
% % % 
% % % for k = 1:numel(beta)
% % %     
% % %     set(H3,'Xdata',circle2_prime{k}(1,:),'Ydata',circle2_prime{k}(2,:))
% % %     set(H4,'Xdata',P_prime{k}(1,:),'Ydata',P_prime{k}(2,:))
% % %     set(U1,'Xdata',alpha(k),'Ydata',Utotal(k))
% % %     
% % %     for ii = 1:length(P_prime{k}(1,:))
% % %         set(U2(ii),'Xdata',alpha(k),'Ydata',Ui(k,ii))
% % %     end
% % %     
% % %     f=getframe(gcf);
% % %     if mCounter == 1
% % %         [mov(:,:,1,mCounter), map] = rgb2ind(f.cdata, 256, 'nodither');
% % %     else
% % %         mov(:,:,1,mCounter) = rgb2ind(f.cdata, map, 'nodither');
% % %     end
% % %     mCounter = mCounter+1;
% % % end
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%               Save animations                  %%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % directory = 'C:\Users\Reinier\Box Sync\Reinier PhD\3_Matlab\SB_Tusi_couple\Animations';
% % % directory = 'D:\pkuppens\Box Sync\Reinier PhD\3_Matlab\SB_Tusi_couple\Animations';
% % % 
% % % idx=0;
% % % cd(directory)
% % % existingNames=dir('*.gif');
% % % for k = 1:numel(existingNames)
% % %     idx(k)=str2num(existingNames(k).name(end-7:end-4));
% % % end
% % % maxIdx = max(idx);
% % % 
% % % animation_name = strcat(directory,'\Animation',num2str(maxIdx+1.','%04d'),'.gif');
% % % imwrite(mov, map,animation_name, 'DelayTime', 0, 'LoopCount', inf);
% % % 
% % % 
% % % 
% % % 
% % % 










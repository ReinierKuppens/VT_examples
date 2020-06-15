
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           Parameter Definitions                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha   = linspace(0,2*pi,200);

r1      = 10;
r2      = 5;
r3      = 4;

nSpring = 3;

rho     =  (r1/r2)-1;
beta    = alpha.*rho;

circle1 = [r1*cos(alpha)        ;r1*sin(alpha)];
circle2 = [r2*cos(alpha)        ;r2*sin(alpha)];
circle3 = [(r1-r2)*cos(alpha)   ;(r1-r2)*sin(alpha)];

theta  = linspace(0,2*pi,nSpring+1);theta(end)=[];
P      = r3.*[cos(theta);sin(theta)];

R  = @(a) [cos(a),-sin(a);sin(a),cos(a)];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%     Compute  distance^2 and energy             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

circle2_prime   = cell(1,numel(alpha));
P_prime         = cell(1,numel(alpha));
%Compute point locations of complete domain alpha
for k = 1:numel(alpha)
    circle2_prime{k}    = R(-alpha(k)*rho)*circle2 + [(r1-r2)*cos(alpha(k));(r1-r2)*sin(alpha(k))];
    P_prime{k}          = R(-alpha(k)*rho)*P       + [(r1-r2)*cos(alpha(k));(r1-r2)*sin(alpha(k))];
end

%Compute energies
stifness = ones(1,nSpring) ;
for k = 1:numel(beta)
    for ii = 1:nSpring
        Ui(k,ii) = 0.5*stifness(ii)*(P_prime{k}(1,ii)^2 + P_prime{k}(2,ii)^2);
    end
end
Utotal = sum(Ui,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                Initial Plots                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
H4 = plot(P(1,:)+circle3(1,1),P(2,:)+circle3(2,1),'b.','markersize',20);

plot([0,0],[-r1,r1],'color',[0.7,0.7,0.7])
plot(circle3(1,:),circle3(2,:),'color',[0.7,0.7,0.7])
xlabel('Distance [m]')
ylabel('Distance [m]')

for k = 1:numel(beta)
    plot(P_prime{k}(1,:),P_prime{k}(2,:),'.r','markersize',5)
end

subplot(1,2,2)

plot(alpha,Utotal);hold on;
for k = 1:length(P_prime{1}(1,:))
    plot(alpha,Ui(:,k));
end

U1 = plot(alpha(1),Utotal(1),'.b','markersize',20);
for ii = 1:length(P_prime{k}(1,:))
    U2(ii) = plot(alpha(1),Ui(1,ii),'.r','markersize',20);
end
xlabel('angle [rad]')
ylabel('Energy [J]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                Make animation                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mCounter=1;
for k = 1:numel(beta)
    
    set(H3,'Xdata',circle2_prime{k}(1,:),'Ydata',circle2_prime{k}(2,:))
    set(H4,'Xdata',P_prime{k}(1,:),'Ydata',P_prime{k}(2,:))
    set(U1,'Xdata',alpha(k),'Ydata',Utotal(k))
    
for ii = 1:length(P_prime{k}(1,:))
    set(U2(ii),'Xdata',alpha(k),'Ydata',Ui(k,ii))
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
















clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           Parameter Definitions                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha   = linspace(0,2*pi,200);

r1      = 10;
r2      = 5;
r3      = 5;

nSpring = 3;

rho     = 1;
beta    = alpha.*rho;

circle1 = [r1*cos(alpha)        ;r1*sin(alpha)];
circle2 = [r2*cos(alpha)        ;r2*sin(alpha)];
circle3 = [(r1-r2)*cos(alpha)   ;(r1-r2)*sin(alpha)];

T = @(a,tx)[cos(a),-sin(a),tx;
            sin(a), cos(a),0 ;
            0     , 0     ,0 ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%     Random         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m1      = 1;
g       = 9.81;

x1 = -r3 + (2*r3).*rand([1,nSpring]);
y1 = -r3 + (2*r3).*rand([1,nSpring]);
r  = r3.*rand([1,nSpring]);

A1          = sqrt(x1.^2 + y1.^2); 
theta       = atan2(y1,x1);
stiffness   = A1./r;

x2 = sum(x1);
y2 = sum(y1);

A2      = sqrt(x2.^2 + y2.^2); 
eta1    = atan2(y2,x2);
q1      = ((r2-r1).*A2)./(g*m1);

check1 = sum(stiffness.*r.*cos(theta)) - g*m1*q1*cos(eta1)./(r2-r1);
check2 = sum(stiffness.*r.*sin(theta)) - g*m1*q1*sin(eta1)./(r2-r1);

% keyboard

P   =   r.*[cos(theta);sin(theta)];
M1  =  q1.*[cos(eta1) ;sin(eta1)];


[x1,y1,x2,y2,x3,y3]=getComOutline([M1(1);M1(2)],0.5);
com1 = [x1;y1];com2=[x2;y2];com3=[x3;y3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%     Compute  distance^2 and energy             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

circle2_prime   = cell(1,numel(alpha));
P_prime         = cell(1,numel(alpha));
M1_prime        = cell(1,numel(alpha));
com1_prime            = cell(1,numel(alpha));
com2_prime            = cell(1,numel(alpha));
com3_prime            = cell(1,numel(alpha));

%Compute point locations of complete domain alpha
for k = 1:numel(alpha)
    
    circle2_prime{k}    = T(alpha(k),0)*T(alpha(k)*rho,(r1-r2))*[circle2;ones(1,length(circle2(1,:)))];
    P_prime{k}          = T(alpha(k),0)*T(alpha(k)*rho,(r1-r2))*[P;ones(1,length(P(1,:)))];

    M1_prime{k}         = T(alpha(k),0)*[M1;1];
    com1_prime{k}       = T(alpha(k),0)*[com1;ones(1,length(com1))];
    com2_prime{k}       = T(alpha(k),0)*[com2;ones(1,length(com2))];
    com3_prime{k}       = T(alpha(k),0)*[com3;ones(1,length(com3))];

 
end

%Compute energies
% stifness = ones(1,nSpring) ;
for k = 1:numel(beta)
    for ii = 1:nSpring
        Ui(k,ii) = 0.5*stiffness(ii)*(P_prime{k}(1,ii)^2 + P_prime{k}(2,ii)^2);
    end
    Ug(k) = m1*g*M1_prime{k}(1);
end
Utotal = sum(Ui,2) + Ug.';
% keyboard 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                Initial Plots                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fh      = figure('color',[1,1,1],'position',[100,100,1000,500]);subplot(1,2,1);hold on ;axis equal;
ht      = title('Gravity pointed in negative x direction'); %Print time in title

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




H5 = fill(com3_prime{k}(1,:),com3_prime{k}(2,:),'w');
H6 = plot(com3_prime{k}(1,:),com3_prime{k}(2,:),'k');
H7 = fill(com1_prime{k}(1,:),com2_prime{k}(2,:),'k');
H8 = fill(com2_prime{k}(1,:),com1_prime{k}(2,:),'k');

% H5 = plot(M1_prime{1}(1),M1_prime{1}(2),'k.','markersize',10)


plot([0,0],[-r1,r1],'color',[0.7,0.7,0.7])
plot(circle3(1,:),circle3(2,:),'color',[0.7,0.7,0.7])
xlabel('Distance [m]')
ylabel('Distance [m]')

for k = 1:numel(beta)
    plot(P_prime{k}(1,:),P_prime{k}(2,:),'.k','markersize',5)
end

subplot(1,2,2)

plot(alpha,Utotal,'r','linewidth',2);hold on;
for k = 1:length(P_prime{1}(1,:))
    plot(alpha,Ui(:,k));
end
plot(alpha,Ug,'--k','linewidth',2)

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
    
    
    [x1,y1,x2,y2,x3,y3]=getComOutline([M1_prime{k}(1);M1_prime{k}(2)],1);

    set(H5,'Xdata',com3_prime{k}(1,:),'Ydata',com3_prime{k}(2,:))
    set(H6,'Xdata',com3_prime{k}(1,:),'Ydata',com3_prime{k}(2,:))
    set(H7,'Xdata',com1_prime{k}(1,:),'Ydata',com1_prime{k}(2,:))
    set(H8,'Xdata',com2_prime{k}(1,:),'Ydata',com2_prime{k}(2,:))

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















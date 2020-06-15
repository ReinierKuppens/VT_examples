

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           Parameter Definitions                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r1 = 10;
r2 = 5;

nSpring = 1;

rho1 = 1;
rho2 = 1;
rho3 = 1;

alpha  = linspace(pi,3*pi,200);
beta   = linspace(0,3*pi,200);
gamma  = linspace(pi,3*pi,200);
 
[X1,Y1,Z1] = sphere(50);
[X2,Y2,Z2] = sphere(50);
[X3,Y3,Z3] = sphere(50);

sx1 = X1*r1;
sy1 = Y1*r1;
sz1 = Z1*r1;
sx2 = X2*r2;
sy2 = Y2*r2;
sz2 = Z2*r2;
sx3 = X3*(r1-r2);
sy3 = Y3*(r1-r2);
sz3 = Z3*(r1-r2);
sArrayShape = size(sx2);

sphere2 = [sx2(:).' ; sy2(:).' ; sz2(:).' ; ones(1,numel(sx2))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           Spring and mass postion              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


r3 = 10;

x3 = r3*rand(1);


x1 = -r3 + (2*r3).*rand([1,nSpring]);
y1 = -r3 + (2*r3).*rand([1,nSpring]);
z1 = -r3 + (2*r3).*rand([1,nSpring]);
r  =  r3.*rand([1,nSpring]);

A1          = sqrt(x1.^2 + y1.^2 + z1.^2);
theta       = acos(z1./A1);
phi         = atan2(y1,x1);

stiffness   = A1./r;

g   = 9.81;
m1  = 1;

x2 = sum(x1)+x3;
y2 = sum(y1);
z2 = sum(z1);

A2       = sqrt(x2.^2 + y2.^2 + z2.^2);
eta1     = acos(z2./A2);
zeta1    = atan2(y2,x2);

q1      = ((r2-r1).*A2)./(g.*m1);

m2      = x3./g; 
q2      = 0; 
eta2    = 0;
zeta2   = 0; 

check1 = sum(stiffness.*r.*(r2-r1).*sin(theta).*cos(phi))   - g.*m1.*q1.*sin(eta1).*cos(zeta1) + g.*m2.*(r2-r1);
check2 = sum(stiffness.*r.*(r2-r1).*sin(theta).*sin(phi))   - g.*m1.*q1.*sin(eta1).*sin(zeta1);
check3 = sum(stiffness.*r.*(r2-r1).*cos(theta))             - g.*m1.*q1.*cos(eta1);


P = [   r.*sin(theta).*cos(phi)
        r.*sin(theta).*sin(phi)
        r.*cos(theta)
        ones(size(theta))];


M1 = [  q1.*sin(eta1).*cos(zeta1)
        q1.*sin(eta1).*sin(zeta1)
        q1.*cos(eta1)
        ones(size(eta1))];
    
M2 = [  q2.*sin(eta2).*cos(zeta2)
        q2.*sin(eta2).*sin(zeta2)
        q2.*cos(eta2)
        ones(size(eta1))];
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%       def Homogeneous tranformations           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rx = @(a) [ 1       0       0 ;
            0       cos(a)  -sin(a);
            0       sin(a)  cos(a);];

Ry = @(a) [ cos(a)  0       sin(a);
            0       1       0;
           -sin(a)  0       cos(a)];
        
Rz = @(a) [ cos(a) -sin(a)  0;
            sin(a)  cos(a)  0;
            0       0       1];


Rzyx = @(a,b,c) Rx(a)*Ry(b)*Rz(c);
H    = @(a,b,c,t) [Rzyx(a,b,c), [t;0;0] ; [0,0,0,1]];



[Xcom1,Ycom1,Zcom1,Ccom1]=getComOutline3D(M1,0.5); 
[Xcom2,Ycom2,Zcom2,Ccom2]=getComOutline3D(M2,0.5); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%        Compute rotations and energies          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M1_prime    = cell(1,numel(alpha));
M2_prime    = cell(1,numel(alpha));

P_prime     = cell(1,numel(alpha));
sphere2_prime = cell(1,numel(alpha));
Com1_prime = cell(1,numel(alpha));
Com2_prime = cell(1,numel(alpha));

for k = 1:numel(gamma)
    
    a1 = alpha(k);
    b1 = beta(k);
    c1 = gamma(k);
    t1 = 0;
    
    a2 = rho1.*alpha(k);
    b2 = rho2.*beta(k);
    c2 = rho3.*gamma(k);
    t2 = r1-r2;
    
    P_prime{k}          = H(a1,b1,c1,t1)*H(a2,b2,c2,t2)*P;
    sphere2_prime{k}    = H(a1,b1,c1,t1)*H(a2,b2,c2,t2)*sphere2;
    
    M1_prime{k}         = H(a1,b1,c1,t1)*M1;
    M2_prime{k}         = H(a1,b1,c1,t1)*H(a2,b2,c2,t2)*M2;

    Com1_prime{k}        = H(a1,b1,c1,t1)*[Xcom1(:).';Ycom1(:).';Zcom1(:).';ones(1,numel(Xcom1(:)))];
    Com2_prime{k}        = H(a1,b1,c1,t1)*H(a2,b2,c2,t2)*[Xcom2(:).';Ycom2(:).';Zcom2(:).';ones(1,numel(Xcom2(:)))];
end

for k = 1:numel(gamma)
    for ii = 1:nSpring
        Ui(k,ii) = 0.5*stiffness(ii)*(P_prime{k}(1,ii)^2 + P_prime{k}(2,ii)^2 + P_prime{k}(3,ii)^2);
    end
    Ug1(k) = m1*g*M1_prime{k}(1);
    Ug2(k) = m2*g*M2_prime{k}(1);
end

UtotalS = sum(Ui,2); 
Utotal  = sum(Ui,2) + Ug1.' + Ug2.';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                Initial Plots                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh      = figure('color',[1,1,1],'position',[100,100,1000,500]);subplot(1,2,1);hold on ;axis equal;
ht      = title(['R1 = ',num2str(r1),', R2 = ',num2str(r2),', nSpring = ',num2str(nSpring)]); %Print time in title
xlabel('x-axis');ylabel('y-axis');zlabel('z-axis')

for k = 1:nSpring
    for ii = 1:numel(P_prime)
        Px(k,ii) = P_prime{ii}(1,k);
        Py(k,ii) = P_prime{ii}(2,k);
        Pz(k,ii) = P_prime{ii}(3,k);
    end
end

for k = 1:numel(M1_prime)
    M1x(k) = M1_prime{k}(1);
    M1y(k) = M1_prime{k}(2);
    M1z(k) = M1_prime{k}(3);
    M2x(k) = M2_prime{k}(1);
    M2y(k) = M2_prime{k}(2);
    M2z(k) = M2_prime{k}(3);
end


% keyboard

f       = getframe(gcf);
sizef   = size(f.cdata);
width   = sizef(1);
height  = sizef(2);


subplot(1,2,1); hold on
material shiny

surface(sx1,sy1,sz1,'EdgeAlpha', 0.1,...
    'FaceAlpha', 0.3,...
    'FaceColor', [189, 189, 189]./255,...
    'EdgeColor', [189, 189, 189]./255,...
    'AmbientStrength',0.5)

surface(sx3,sy3,sz3,'EdgeAlpha', 0.2,...
    'FaceAlpha', 0.8,...
    'FaceColor', [0 188  212]./255,...
    'EdgeColor', 'none')

sx2_prime = reshape(sphere2_prime{1}(1,:),sArrayShape(1),sArrayShape(2));
sy2_prime = reshape(sphere2_prime{1}(2,:),sArrayShape(1),sArrayShape(2));
sz2_prime = reshape(sphere2_prime{1}(3,:),sArrayShape(1),sArrayShape(2));

H1 = surface(sx2_prime,sy2_prime,sz2_prime,...
    'FaceAlpha', 0.4, ...
    'FaceColor', [0.1 0.1 0.1],...
    'EdgeColor', 'none',...
    'FaceLighting','gouraud',...
    'AmbientStrength',0.5);

H2 = plot3(P_prime{1}(1,:),P_prime{1}(2,:),P_prime{1}(3,:),'.k','markersize',20);

H5 = surface(Xcom1,Ycom1,Zcom1,'Cdata',Ccom1,...
    'FaceAlpha', 1, ...
    'FaceLighting','gouraud',...
    'AmbientStrength',0.5,...
    'EdgeColor', 'none');

H6 = surface(Xcom2,Ycom2,Zcom2,'Cdata',Ccom2,...
    'FaceAlpha', 1, ...
    'FaceLighting','gouraud',...
    'AmbientStrength',0.5,...
    'EdgeColor', 'none');


dist = 12;
patch([dist -dist -dist dist], [dist  dist -dist -dist],  [-dist -dist -dist -dist], [120, 144, 156]./255,'facealpha',0.5)
patch([dist -dist -dist dist], [dist  dist  dist  dist],  [ dist  dist -dist -dist], [120, 144, 156]./255,'facealpha',0.5)
patch([dist  dist  dist dist], [dist -dist -dist  dist],  [ dist  dist -dist -dist], [120, 144, 156]./255,'facealpha',0.5)

light('Position',[10 -10 10],'Style','local')
% alpha 0.5

for k = 1:nSpring
    plot3(Px(k,:),Py(k,:),Pz(k,:),'-r','Linewidth',2)
end
plot3(M1x,M1y,M1z,'--k')
plot3(M2x,M2y,M2z,':k')


axis equal
view(3)

subplot(1,2,2);hold on
U1 = plot(gamma(:),Utotal(:),'-b');
for ii = 1:length(P_prime{k}(1,:))
    U2(ii) = plot(gamma(:),Ui(:,ii),'-r');
end
plot(gamma(:),Ug1,'--k')
plot(gamma(:),Ug2,':k')

plot(gamma(:),UtotalS,'--r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                Make animation                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% keyboard 
mCounter=1;

for k = 1:numel(gamma)
    
    sx2_prime = reshape(sphere2_prime{k}(1,:),sArrayShape(1),sArrayShape(2));
    sy2_prime = reshape(sphere2_prime{k}(2,:),sArrayShape(1),sArrayShape(2));
    sz2_prime = reshape(sphere2_prime{k}(3,:),sArrayShape(1),sArrayShape(2));
    
    
%     [Xcom,Ycom,Zcom,Ccom]=getComOutline3D(M1_prime{k},1);

Xcom1_prime = reshape(Com1_prime{k}(1,:),size(Xcom1));
Ycom1_prime = reshape(Com1_prime{k}(2,:),size(Xcom1));
Zcom1_prime = reshape(Com1_prime{k}(3,:),size(Xcom1));

Xcom2_prime = reshape(Com2_prime{k}(1,:),size(Xcom2));
Ycom2_prime = reshape(Com2_prime{k}(2,:),size(Xcom2));
Zcom2_prime = reshape(Com2_prime{k}(3,:),size(Xcom2));

    
    set(H1, 'Xdata',sx2_prime,...
            'Ydata',sy2_prime,...
            'Zdata',sz2_prime)
    
    
    set(H2, 'Xdata',P_prime{k}(1,:),...
            'Ydata',P_prime{k}(2,:),...
            'Zdata',P_prime{k}(3,:))
    
%     set(H3, 'Xdata',mean(sx2_prime(:)),...
%             'Ydata',mean(sy2_prime(:)),...
%             'Zdata',mean(sz2_prime(:)))
        
%     set(H4, 'Xdata',M2_prime{k}(1),...
%             'Ydata',M2_prime{k}(2),...
%             'Zdata',M2_prime{k}(3))
    
    set(H5, 'Xdata',Xcom1_prime,...
            'Ydata',Ycom1_prime,...
            'Zdata',Zcom1_prime,...            
            'Cdata',Ccom1)
        
    set(H6, 'Xdata',Xcom2_prime,...
            'Ydata',Ycom2_prime,...
            'Zdata',Zcom2_prime,...            
            'Cdata',Ccom1)

    
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
















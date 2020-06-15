

clear all 
close all 
clc 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           Parameter Definitions                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r1 = 10;
r2 = 4;
r3 = 4;

nSpring = 4; 

rho1    = 2;
rho2    = 3;
rho3    = 0.5; 

alpha   = linspace(0,2*pi,200);
gamma   = linspace(0,2*pi,200);
beta    = linspace(0,3*pi,200); 

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

theta = linspace(0,pi,nSpring+1);theta(end)=[];
phi   = linspace(0 ,2*pi,nSpring+1);phi(end)=[];

P = [   r3*sin(theta).*cos(phi)
        r3*sin(theta).*sin(phi)
        r3*cos(theta)
        ones(size(theta))];
    
P = r3.*rand(3,nSpring-1);
P(:,end+1) = -sum(P,2);
P = [P;ones(size(theta))];
    
    
 Rz = @(a) [cos(a) -sin(a)  0; 
            sin(a)  cos(a)  0; 
            0       0       1];
        
 Ry = @(a) [cos(a)  0       sin(a); 
            0       1       0; 
            -sin(a) 0       cos(a)];
        
 Rx = @(a) [1       0       0 ; 
            0       cos(a)  -sin(a); 
            0       sin(a)  cos(a);];
        

 
 Rzyx = @(a,b,c) Rx(a)*Ry(b)*Rz(c)
 
 H = @(a,b,c,t) [Rzyx(a,b,c), [t;0;0] ; [0,0,0,1]]
 
 for k = 1:numel(gamma)
    
    a1 = alpha(k);
    b1 = beta(k);
    c1 = gamma(k); 
    t1 = 0;   
    
    a2 = rho1.*alpha(k);
    b2 = rho2.*beta(k);
    c2 = rho3.*gamma(k);
    t2 = r1-r2;  
    
    P_prime{k} = H(a1,b1,c1,t1)*H(a2,b2,c2,t2)*P;
    sphere2_prime{k} = H(a1,b1,c1,t1)*H(a2,b2,c2,t2)*sphere2;
     
 end
 
 stifness = ones(1,nSpring) ;
for k = 1:numel(gamma)
    for ii = 1:nSpring
        Ui(k,ii) = 0.5*stifness(ii)*(P_prime{k}(1,ii)^2 + P_prime{k}(2,ii)^2 + P_prime{k}(3,ii)^2);
    end
end
 
  Utotal = sum(Ui,2);
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                Initial Plots                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
fh      = figure('color',[1,1,1],'position',[100,100,1000,500]);subplot(1,2,1);hold on ;axis equal;
ht      = title(['R1 = ',num2str(r1),', R2 = ',num2str(r2),', nSpring = ',num2str(nSpring)]); %Print time in title

for k = 1:nSpring
    for ii = 1:numel(P_prime) 
        
        Px(k,ii) = P_prime{ii}(1,k);
        Py(k,ii) = P_prime{ii}(2,k);
        Pz(k,ii) = P_prime{ii}(3,k);

        
    end
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
H3 = plot3(mean(sx2_prime(:)),mean(sy2_prime(:)),mean(sz2_prime(:)),'.g','markersize',10);

dist = 11

patch([dist -dist -dist dist], [dist dist -dist -dist], [-dist -dist -dist -dist],[120, 144, 156]./255,'facealpha',1)
patch([dist -dist -dist dist], [dist dist dist dist], [dist dist -dist -dist], [120, 144, 156]./255,'facealpha',1)
patch([dist dist dist dist],[dist -dist -dist dist],  [dist dist -dist -dist], [120, 144, 156]./255,'facealpha',1)

                    light('Position',[10 -10 10],'Style','local')
%                     shading flat
                    
%                     lighting gouraud


for k = 1:nSpring
   plot3(Px(k,:),Py(k,:),Pz(k,:),'-r','Linewidth',2) 
end

% alpha 0.1 

axis equal
view(3)

subplot(1,2,2);hold on 
  U1 = plot(gamma(:),Utotal(:),'-b');
for ii = 1:length(P_prime{k}(1,:))
    U2(ii) = plot(gamma(:),Ui(:,ii),'-r');
end
       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                Make animation                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mCounter=1;

for k = 1:numel(gamma)

    sx2_prime = reshape(sphere2_prime{k}(1,:),sArrayShape(1),sArrayShape(2));
    sy2_prime = reshape(sphere2_prime{k}(2,:),sArrayShape(1),sArrayShape(2));
    sz2_prime = reshape(sphere2_prime{k}(3,:),sArrayShape(1),sArrayShape(2));
    
    set(H1, 'Xdata',sx2_prime,...
            'Ydata',sy2_prime,...
            'Zdata',sz2_prime)
    
    
    set(H2, 'Xdata',P_prime{k}(1,:),...
            'Ydata',P_prime{k}(2,:),...
            'Zdata',P_prime{k}(3,:))
        
    set(H3, 'Xdata',mean(sx2_prime(:)),...
            'Ydata',mean(sy2_prime(:)),...
            'Zdata',mean(sz2_prime(:)))
    
%         pause(0.05)
        
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
















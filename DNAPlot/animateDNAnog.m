%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %
%   Original Date: 29-06-2015                                             %
%
%   Last changed: 2019-08-24
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%animateDNA(DNA,varargin) generates an animation of the mechanism described by DNA.


function animateDNAnog(DNA,t,qss)

fh = figure('color',[1,1,1],'position',[100,100,1000,500]);subplot(1,2,1);hold on ;axis equal;

[massPoints,comPoint] = getMassLocations(DNA);
nM                    = numel(massPoints);
R                     = @(angle) [cos(angle),-sin(angle);sin(angle),cos(angle)];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%                Initial Plots                   %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%                Plot Grounds                   %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub1h = subplot(1,2,1);
sub1h.Position = sub1h.Position + [-0.1 -0.05 0.1 0.1];

hold on
axis([-3.2 3.2 -3.2 3.2])
xlabel('Distance [m]')
ylabel('Distance [m]')
axis off

% One arrow from left to right with text on left side
if ~any(DNA.Mpar(:,3)==0)
    xTEXT = [0.1 0.05];       %length and location of arrow
    yTEXT = [0.85 0.85];      %height and width of arrow
    annotation('textarrow',xTEXT,yTEXT,'String','g','FontSize',13,'Linewidth',2,'interpreter','latex')
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%                Plot Trajectories               %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[allTrajEquallySpaced,nTraj] = getAllTrajectories(DNA,t,qss);
for k = 6:nTraj
    plot(allTrajEquallySpaced(:,k*2-1),allTrajEquallySpaced(:,k*2),'k')
end


plotGrounds(DNA);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%                Plot Bodies                     %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cdata   = rand(nM,3);
% cdata   = [209 237 243;209 237 243]./255;

%Get the initial outline of the parts with the centre of mass point

for k = 1:nM
    [x1{k},y1{k},x2{k},y2{k},x3{k},y3{k}]   = getComOutline(comPoint{k},0.075); %ROW VECTORS
    [polx{k},poly{k},scx{k},scy{k}]         = getMassOutline(massPoints{k}); %ROW VECTORS
end

%Draw all masses and create plot handles
for k =nM:-1:1
    polh(k)=fill(polx{k},poly{k},cdata(k,:),'Edgecolor','k');
    alpha(0.5)
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%                Plot COMS                      %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

comColors = {'k','b','r','m','c'};


%Draw all com points and create plot handles

for k =1:nM
    for l =1:numel(scx{k})
        polh2(k,l)=fill(scx{k}{l},scy{k}{l},'w','Linewidth',1);
    end
        
        comh1(k)=fill(x3{k},y3{k},'w');
        comh2(k)=plot(x3{k},y3{k},comColors{k});
        comh3(k)=fill(x1{k},y1{k},comColors{k});
        comh4(k)=fill(x2{k},y2{k},comColors{k});
    end
    





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%                Plot Springs                    %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



sSDNA   = size(DNA.Spar);
nS      = sSDNA(2);
incS    = DNA2inc(DNA,'S');

if nS>0
    EsColor = getSpringColor(DNA,t,qss);
end

if nS>0
    [xs,ys,nCoils,maxL,distUncoiled] = getSpringOutline(DNA,t,qss);
    for k = 1:nS
        cdataS = EsColor{k}(1,:);
        springh(k) = plot(xs{k},ys{k},'color',cdataS,'Linewidth',1);
    end
end

axis square
axis equal

% keyboard

ylim(sub1h,sub1h.YLim+[-0.1,0.1])  % restore y limits
xlim(sub1h,sub1h.XLim+[-0.1,0.1])

% keyboard
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%                Plot Energies                  %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sub2h = subplot(1,2,2);
% sub2h.Position = sub2h.Position + [0 0 -0.03 0];

hold on
xlabel('\alpha [rad]')
xticks([0,0.5*pi,1*pi,1.5*pi,2*pi])
xticklabels({'0','0.5\pi','\pi','1.5\pi','2\pi'})
ylabel('Energy [Nm]')
grid on

[Es,~]  = getEnergies(DNA,t,qss);
Eg      = getGravityEnergy(DNA,t,qss);
% keyboard 

Et = zeros(1,numel(t));
for k = 1:numel(Es)
    energyhAll(k) = patch([t,NaN],[Es{k},NaN],[Es{k},NaN]);
    
    Et = Et + Es{k};
end

% keyboard
Et = Et + sum(Eg,1);

for k = 1:numel(Es)
    energyh(k)    = plot(t(1),Es{k}(1),'.r','markersize',30);
    
    lineColor = [EsColor{k};EsColor{k}(end,:)];
    set(energyhAll(k),'edgecolor','interp','linewidth',2,'FaceVertexCData',lineColor)
end


for k = 1:nM
    plot(t,Eg(k,:),'--','linewidth',2,'color',comColors{k});
    gravityh(k)    = plot(t(1),Eg(k,1),'.','markersize',30,'color',comColors{k});
    
end

plot(t,Et,'k','Linewidth',2)
totalEnergyh = plot(t(1),Et(1),'.k','markersize',30);
%
% for k = 1:numel(Es)
%     plot(t,Es{k},'k')
% end

% keyboard
axis square
yl=ylim(sub2h); % retrieve auto y-limits
axis tight   % set tight range
ylim(sub2h,yl)  % restore y limits

% keyboard




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%               Start Animation                  %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% legend('Total energy','Individual spring energy','location','northoutside')

mCounter  = 1;

for ii = 1:length(t)
    
    %     set(ht, 'String', ['\rho = ',num2str(-DNA.rho),', \alpha = ',sprintf('%0.2f',t(ii)/(pi)),'\pi rad'])
    %Move the masses
    
    for k =1:nM
        RpolhData  = R(qss(ii,k*3))*[polx{k}-comPoint{k}(1);poly{k}-comPoint{k}(2)];
        RTpolhData = bsxfun(@plus,RpolhData,[qss(ii,k*3-2);qss(ii,k*3-1)]);
        set(polh(k),'Xdata',RTpolhData(1,:),'Ydata',RTpolhData(2,:))
    end
    
    %Move the grounds
    
    for k =1:nM
        Rx1Data = R(qss(ii,k*3))*[x1{k}-comPoint{k}(1);y1{k}-comPoint{k}(2)];
        RTx1Data = bsxfun(@plus,Rx1Data,[qss(ii,k*3-2);qss(ii,k*3-1)]);
        
        Rx2Data = R(qss(ii,k*3))*[x2{k}-comPoint{k}(1);y2{k}-comPoint{k}(2)];
        RTx2Data = bsxfun(@plus,Rx2Data,[qss(ii,k*3-2);qss(ii,k*3-1)]);
        
        Rx3Data = R(qss(ii,k*3))*[x3{k}-comPoint{k}(1);y3{k}-comPoint{k}(2)];
        RTx3Data = bsxfun(@plus,Rx3Data,[qss(ii,k*3-2);qss(ii,k*3-1)]);
        
        
            
            set(comh1(k),'Xdata',RTx3Data(1,:),'Ydata',RTx3Data(2,:))
            set(comh2(k),'Xdata',RTx3Data(1,:),'Ydata',RTx3Data(2,:))
            set(comh3(k),'Xdata',RTx1Data(1,:),'Ydata',RTx1Data(2,:))
            set(comh4(k),'Xdata',RTx2Data(1,:),'Ydata',RTx2Data(2,:))
        
        for l =1:numel(scx{k})
            Rpolh2Data{l}= R(qss(ii,k*3))*[scx{k}{l}-comPoint{k}(1);scy{k}{l}-comPoint{k}(2)];
            RTpolh2Data{l} = bsxfun(@plus,Rpolh2Data{l},[qss(ii,k*3-2);qss(ii,k*3-1)]);
            
            set(polh2(k,l),'Xdata',RTpolh2Data{l}(1,:),'Ydata',RTpolh2Data{l}(2,:))
        end
    end
    
    for k=1:nS
        
        cdataS = EsColor{k}(ii,:);
        massnr = find(incS(:,k)==1);    %Find entries in column
        
        if any(massnr==1)
            mnr = massnr(2)-1;
            
            com            = comPoint{mnr}.';
            origo{k}       = DNA.Spar(1:2,k);
            insertie{k}    = R(qss(ii,mnr*3))*(DNA.Spar(3:4,k)-com) + [qss(ii,mnr*3-2);qss(ii,mnr*3-1)];
            
            [x{k} y{k}] = springcoord(origo{k},insertie{k},nCoils(k),maxL(k),distUncoiled(k));
            
            set(springh(k),'Xdata',x{k},'Ydata',y{k},'Color',cdataS)
            
        else
            mnr1 = massnr(1)-1;
            mnr2 = massnr(2)-1;
            
            com1 = comPoint{mnr1}.';
            com2 = comPoint{mnr2}.';
            
            origo{k}    = R(qss(ii,mnr1*3))*(DNA.Spar(1:2,k)-com1)+[qss(ii,mnr1*3-2);qss(ii,mnr1*3-1)];
            insertie{k} = R(qss(ii,mnr2*3))*(DNA.Spar(3:4,k)-com2)+[qss(ii,mnr2*3-2);qss(ii,mnr2*3-1)];
            
            [x{k} y{k}] = springcoord(origo{k},insertie{k},nCoils(k),maxL(k),distUncoiled(k));
            
            set(springh(k),'Xdata',x{k},'Ydata',y{k},'Color',cdataS)
            
        end
        
    end
    
    for k = 1:numel(Es)
        cdataS = EsColor{k}(ii,:);
        
        set(energyh(k),'Xdata',t(ii),'Ydata',Es{k}(ii),'Color',cdataS)
    end
    set(totalEnergyh,'Xdata',t(ii),'Ydata',Et(ii))
    
%     if ~any(DNA.Mpar(:,3)==0)

    for k = 1:nM
        
        set(gravityh(k),'Xdata',t(ii),'Ydata',Eg(k,ii),'Color',comColors{k})
    end
%     end
    %     keyboard
    
    f=getframe(gcf);
    
    if mCounter == 1
        [mov(:,:,1,mCounter), map] = rgb2ind(f.cdata, 256, 'nodither');
    else
        mov(:,:,1,mCounter) = rgb2ind(f.cdata, map, 'nodither');
    end
    
    imwrite(f.cdata,['C:\Users\pkuppens\Box Sync\Reinier PhD\3_Matlab\SB_Tusi_couple\Animations\SeparateFiles\ani',num2str(mCounter),'.png'])
    
    mCounter = mCounter+1;
    
    
    
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%              Save Animation to disk            %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

directory = 'C:\Users\pkuppens\Box Sync\Reinier PhD\3_Matlab\SB_Tusi_couple\Animations';

idx=0;
cd(directory)
existingNames=dir('*.gif');
for k = 1:numel(existingNames)
    idx(k)=str2num(existingNames(k).name(end-7:end-4));
end
maxIdx = max(idx);

% keyboard

animation_name = strcat(directory,'\Animation',num2str(maxIdx+1.','%04d'),'.gif');
imwrite(mov, map,animation_name, 'DelayTime', 0, 'LoopCount', inf);






end%Function





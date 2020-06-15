%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %
%   Original Date: 29-06-2015                                             %
%
%   Last changed: 2019-08-24
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%animateDNA(DNA,varargin) generates an animation of the mechanism described by DNA.


function animateDNA(DNA,t,qss)

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

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
sub1h.Position = sub1h.Position + [-0.1 -0.1 0.1 0.1];

hold on
xlabel('Distance [m]')
ylabel('Distance [m]')
axis off


% axis([0,1,0,1])

% keyboard

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%                Plot Trajectories               %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[allTrajEquallySpaced,nTraj] = getAllTrajectories(DNA,t,qss);
for k = 1:nTraj
    plot(allTrajEquallySpaced(:,k*2-1),allTrajEquallySpaced(:,k*2),'w')
end

comColors = {'--k','--b','--r','--m','--c'};
for k = 1:nM
    plot(qss(:,k*3-2),qss(:,k*3-1),comColors{k},'linewidth',2)
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%             Get max x and y limits             %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

springEndPointsX = DNA.Spar(1:2:4,:);
springEndPointsY = DNA.Spar(2:2:4,:);

xmin = min(springEndPointsX(:));
xmax = max(springEndPointsX(:));
ymin = min(springEndPointsY(:));
ymax = max(springEndPointsY(:));

for k = 1:nTraj
    
    xmin = min([allTrajEquallySpaced(:,k*2-1);xmin]);
    xmax = max([allTrajEquallySpaced(:,k*2-1);xmax]);
    ymin = min([allTrajEquallySpaced(:,k*2);ymin]);
    ymax = max([allTrajEquallySpaced(:,k*2);ymax]);
end

daxis = 0.1;
axis([xmin-daxis, xmax+daxis, ymin-daxis, ymax+daxis])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%           Gravity Arrow                        %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~any(DNA.Mpar(:,3)==0)
    xTEXT = [0.1 0.05];       %length and location of arrow
    yTEXT = [0.85 0.85];      %height and width of arrow
    annotation('textarrow',xTEXT,yTEXT,'String','g','FontSize',20,'Linewidth',2,'interpreter','latex')
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%                PlotGrounds                     %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotGrounds(DNA);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%                Plot Bodies                     %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cdata   = rand(nM,3);
% cdata = linspecer(1);

if nM <= 3
    cdata= [0.1       0.1       0.1;
        0.3467    0.5360    0.6907;
        0.9153    0.2816    0.2878];
else
    cdata = linspecer(nM);
end

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
        springh(k)      = plot(xs{k},ys{k},'color',cdataS,'Linewidth',1);
        
        
        radius2     = 0.02;         %Radius of each hinge axle
        idxcircle   = 0:0.1:2*pi;

        smallcirclex = radius2*cos(idxcircle);
        smallcircley = radius2*sin(idxcircle);
        
        springendh1(k) = fill(smallcirclex+xs{k}(1),smallcircley+ys{k}(1),cdataS,'edgecolor','none')
        springendh2(k) = fill(smallcirclex+xs{k}(end),smallcircley+ys{k}(end),cdataS,'edgecolor','none')
        
        
%         springendh(k)   = plot([xs{k}(1),xs{k}(end)],[ys{k}(1),ys{k}(end)],'.','markersize',20,'color',cdataS);
    end
end

axis square
axis equal

ylim(sub1h,sub1h.YLim+[-0.1,0.1])  % restore y limits
xlim(sub1h,sub1h.XLim+[-0.1,0.1])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%                Plot Energies                  %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Labelfontsize=16;
sub2h = subplot(1,2,2);
% sub2h.Position = sub2h.Position + [0 0 -0.03 0];
% sub2h = figure;
hold on
xlabel('$\alpha$ [rad]','FontSize', Labelfontsize)
xticks([0,0.5*pi,1*pi,1.5*pi,2*pi])
xticklabels({'$0$','$0.5\pi$','$\pi$','$1.5\pi$','$2\pi$'})
ylabel('Energy [Nm]','FontSize', Labelfontsize)
grid on
set(gca,'Fontsize',16,'XColor','k','YColor','k');



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

axis square
yl=ylim(sub2h); % retrieve auto y-limits
axis tight      % set tight range
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
            set(springendh1(k),'Xdata',smallcirclex+x{k}(1),'Ydata',smallcircley+y{k}(1),'FaceColor',cdataS)
            set(springendh2(k),'Xdata',smallcirclex+x{k}(end),'Ydata',smallcircley+y{k}(end),'FaceColor',cdataS)
    end
    
    for k = 1:numel(Es)
        cdataS = EsColor{k}(ii,:);
        
        set(energyh(k),'Xdata',t(ii),'Ydata',Es{k}(ii),'Color',cdataS)
    end
    set(totalEnergyh,'Xdata',t(ii),'Ydata',Et(ii))
    
    
    for k = 1:nM
        
        set(gravityh(k),'Xdata',t(ii),'Ydata',Eg(k,ii),'Color',comColors{k})
    end
    
    
    f=getframe(gcf);
    
%     videoFrame{mCounter} = f.cdata;
    
    
    if mCounter == 1
        [mov(:,:,1,mCounter), map] = rgb2ind(f.cdata, 256, 'nodither');
    else
        mov(:,:,1,mCounter) = rgb2ind(f.cdata, map, 'nodither');
    end
    
    imwrite(f.cdata,['C:\Users\pkuppens\Box Sync\Reinier PhD\3_Matlab\SB_Tusi_couple\Animations\SeparateFiles\ani',num2str(mCounter),'.png'])
    
    print('-dpng','-r300','plot.png');
    videoFrame{mCounter}=imread('plot.png');
    
    mCounter = mCounter+1;
    
    
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%              Save Animation to disk            %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currentDir = cd;

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
video_name     = strcat(directory,'\Animation',num2str(maxIdx+1.','%04d'));


imwrite(mov, map,animation_name, 'DelayTime', 0, 'LoopCount', inf);


% keyboard 
% video = VideoWriter(video_name,'Uncompressed AVI'); %create the video object
video = VideoWriter(video_name,'MPEG-4'); %create the video object


nFirstFrame = 10; 
nRepetitions = 2; 

open(video); %open the file for writing
for k=1:numel(videoFrame) %where N is the number of images
    
    if k == 1
        for ii = 1:nFirstFrame
            writeVideo(video,videoFrame{k}); %write the image to file
        end
    end
    
    for ii = 1:nRepetitions
        writeVideo(video,videoFrame{k}); %write the image to file
    end
end
close(video); %close the file



cd(currentDir);



end%Function





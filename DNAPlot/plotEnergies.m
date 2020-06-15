

function []=plotEnergies(DNA,t,qss)


set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


Labelfontsize = 16; 
sub2h = figure('color',[1,1,1]);
set(gca,'FontSize', Labelfontsize)

comColors = {'k','b','r','m','c'};
[massPoints,comPoint] = getMassLocations(DNA);
nM                    = numel(massPoints);


sSDNA   = size(DNA.Spar);
nS      = sSDNA(2);
incS    = DNA2inc(DNA,'S');

if nS>0
    EsColor = getSpringColor(DNA,t,qss);
end




% xlabel('Displacement [mm]','interpreter','latex','FontSize', Labelfontsize)
% ylabel('Force [N]','interpreter','latex','FontSize', Labelfontsize)


% set(gca,'Fontsize',16,'XColor','k','YColor','k');

hold on
xlabel('$\alpha$ [rad]','FontSize', Labelfontsize)
xticks([0,0.5*pi,1*pi,1.5*pi,2*pi])
xticklabels({'$0$','$0.5\pi$','$\pi$','$1.5\pi$','$2\pi$'})
ylabel('Energy [Nm]','FontSize', Labelfontsize)
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
%     energyh(k)    = plot(t(1),Es{k}(1),'.r','markersize',30);
    
    lineColor = [EsColor{k};EsColor{k}(end,:)];
    set(energyhAll(k),'edgecolor','interp','linewidth',2,'FaceVertexCData',lineColor)
end


for k = 1:nM
    plot(t,Eg(k,:),'--','linewidth',2,'color',comColors{k});
%     gravityh(k)    = plot(t(1),Eg(k,1),'.','markersize',30,'color',comColors{k});
    
end

plot(t,Et,'k','Linewidth',2)
% totalEnergyh = plot(t(1),Et(1),'.k','markersize',30);
% keyboard   
axis square
yl=ylim(sub2h.CurrentAxes); % retrieve auto y-limits
axis tight      % set tight range
ylim(sub2h.CurrentAxes,yl)  % restore y limits









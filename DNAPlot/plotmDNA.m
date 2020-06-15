% plotmDNA(DNA) plot a mechanism. 
%
%   Input:  DNA
%   Output: figure with mechanism of DNA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %           
%   Original Date: 2-07-2015                                              %
%   Last modified: 2-07-2015                                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plotmDNA(DNA,t,qss)



figure('color',[1,1,1]); 


if ~any(DNA.Mpar(:,3)==0)
    xTEXT = [0.15 0.05];       %length and location of arrow
    yTEXT = [0.85 0.85];      %height and width of arrow
    annotation('textarrow',xTEXT,yTEXT,'String','g','FontSize',25,'Linewidth',2,'interpreter','latex')
end

% figure_handle = figure('Color',[1 1 1],'Position', [50,50,1000,700],...
%     'visible','on');        %Set background white
hold on                         %Hold on
axis equal                      %Set axis equal
axis off                        %And off
xlabel('Distance [m]')
ylabel('Distance [m]')


nM = length(DNA.Mpar(:,1))-1;


% comColors = {'','','','--k','--k','--k','--k','--k'};
% hqss = hingeState(DNA,t,qss);
% for k = 4:length(hqss(1,:))/2-1
%     plot(hqss(:,k*2-1),hqss(:,k*2),comColors{k},'linewidth',0.5)
% end


comColors = {'--k','--b','--r','--m','--c'};
for k = 1:nM
    plot(qss(:,k*3-2),qss(:,k*3-1),comColors{k},'linewidth',0.5)
end

plotGrounds(DNA);
plotMasses(DNA);
plotSpring(DNA,t,qss);


end%Function











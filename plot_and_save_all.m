
clear all
close all
clc
warning off 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 make names 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name{1}  = '1link_1dof_1spring';
name{2}  = '1link_1dof_3spring';

name{3}  = '2link_1dof_1spring_rho-1_1';
name{4}  = '2link_1dof_1spring_rho-1_2';

name{5}  = '2link_2dof_2spring_1';
name{6}  = '2link_2dof_2spring_2';

name{7}  = '2link_2dof_5spring_1';
name{8}  = '2link_2dof_5spring_2';

name{9}  = '3link_3dof_3spring_1';
name{10} = '3link_3dof_3spring_2';
name{11} = '3link_3dof_3spring_3';
name{12} = '3link_3dof_3spring_4';

name{13}  = '2link_1dof_3spring_rho-1_1';
name{14}  = '2link_1dof_3spring_rho-1_2';

name{15}  = '3link_3dof_2-03_3-12spring';




%loop over all names, load file and save eps and png of mechanism and %energy
for k =  1:15% 13:numel(name)
    close all
    load([name{k},'.mat'])
    
    plotEnergies(DNA,t,state)
    set(gcf, 'Renderer','Painters')
    saveas(gcf,['Pictures/',name{k},'_energy'],'epsc')
    saveas(gcf,['Pictures/',name{k},'_energy'],'png')
    
    
    plotmDNA(DNA,t,state)
    set(gcf, 'Renderer','openGL')
    saveas(gcf,['Pictures/',name{k},'_mech'],'epsc')
    saveas(gcf,['Pictures/',name{k},'_mech'],'png')
    
%     drawnow 
    
    animateDNA(DNA,t,state) 
    
end


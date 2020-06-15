%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %
%   Original Date: 29-06-2015                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [allTrajEquallySpaced,nTraj] = getAllTrajectories(DNA,t,qss)

nM                  = size(DNA.Mpar,1)-1;
qssxy               = qss(:,1:nM*3);        %Duplicate position data
qssxy(:,3:3:nM*3)   = [];                   %Remove angular data

allTraj = [qssxy hingeState(DNA,t,qss)];    %Augment with hinge position data
nTraj   = size(allTraj,2)/2;                %Get number of trajectories



for k = 1:nTraj
    
    if any(abs(sum(diff(allTraj(:,k*2-1:k*2))))>1e-5)
        
        allTrajEquallySpaced(:,k*2-1:k*2) = allTraj(:,k*2-1:k*2);

%         allTrajEquallySpaced(:,k*2-1:k*2) = interparc(lTraj,allTraj(:,k*2-1),allTraj(:,k*2),'lin');
%         allTrajEquallySpaced(:,k*2-1:k*2) = sortPolygonPoints(allTrajEquallySpaced(:,k*2-1:k*2));

    else
        allTrajEquallySpaced(:,k*2-1:k*2) = allTraj(:,k*2-1:k*2);
    end
    
end
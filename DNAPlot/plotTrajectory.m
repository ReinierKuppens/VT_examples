


function [] = plotTrajectory(DNA,t,qss)

figure('color',[1,1,1]);hold on 
axis equal, axis off, 
[allTrajEquallySpaced,nTraj] = getAllTrajectories(DNA,t,qss);


for k = 1:5
    plot(allTrajEquallySpaced(:,k*2-1),allTrajEquallySpaced(:,k*2),'r','Linewidth',1)
end
for k = 6:nTraj
    plot(allTrajEquallySpaced(:,k*2-1),allTrajEquallySpaced(:,k*2),'k','Linewidth',1)
end


% keyboard 
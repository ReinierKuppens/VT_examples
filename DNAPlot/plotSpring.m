


function []=plotSpring(DNA,t,qss)
hold on 


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
        plot(xs{k},ys{k},'color',cdataS,'Linewidth',1);
        
        
        radius2     = 0.025;         %Radius of each hinge axle
        idxcircle   = 0:0.1:2*pi;

        smallcirclex = radius2*cos(idxcircle);
        smallcircley = radius2*sin(idxcircle);
        
        fill(smallcirclex+xs{k}(1),smallcircley+ys{k}(1),cdataS,'edgecolor','none');
        fill(smallcirclex+xs{k}(end),smallcircley+ys{k}(end),cdataS,'edgecolor','none');
        
        
    end
end



% [xs,ys]=getSpringOutline(DNA,t,qss);
% 
% % cdata = 'k';
% cdata = [0 0 0];
% % cdata = [0.5 0.5 0.5];
% % keyboard 
% 
% for k = 1:numel(xs)
%     plot(xs{k},ys{k},'color',cdata,'Linewidth',1)
% end




end
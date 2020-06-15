%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %
%   Original Date: 29-06-2015                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xOutline,yOutline,xSmallcircle,ySmallcircle,xSlide,ySlide] = groundCoords(pointlist,type)

xSlide = [];
ySlide = [];

xTriangle = [-.5 .5 0]+pointlist(1);
yTriangle = [-.5 -.5 0 ]+pointlist(2);

t   = 0:0.2:2*pi;
cx  = 0.1*cos(t) + pointlist(1);
cy  = 0.1*sin(t) + pointlist(2);
xSmallcircle = 0.05*cos(t)+ pointlist(1);
ySmallcircle = 0.05*sin(t)+ pointlist(2);

[polyUnionx,polyUniony]=polybool('union',xTriangle,yTriangle,cx,cy);

switch type
    case 'H'
        
        lw = 0.04; %linewidth of the lines below the ground
        glx = [0 -.1 lw-.1 lw]-0.5;
        gly = [-.5 -.6 -.6 -.5];
        cnt = 0;
        %     for l = 1:11
        nStripes  = 6;
        for l = 1:nStripes
            
            
            xBottomStripe = glx+cnt+pointlist(1);
            yBottomStripe = gly+pointlist(2);
            
            [polyUnionx,polyUniony]=polybool('union',polyUnionx,polyUniony,xBottomStripe,yBottomStripe);
            
            cnt = cnt+0.2-(lw/nStripes);
        end%For
        
    case 'P'
        
        xMax =  5;
        xMin = -5;
        yMin = -0.5;
        
        xSlide = [xMin xMax xMax xMin];
        ySlide = [yMin-0.1 yMin-0.1 yMin-0.2 yMin-0.2];
        
        lw = 0.04; %linewidth of the lines below the ground
        glx = [0 -.1 lw-.1 lw]-0.5;
        gly = [-.5 -.6 -.6 -.5]-0.2;
        cnt = 0;
        %     for l = 1:11
        nStripes  = 6;
        for l = 1:nStripes
            
            
            xBottomStripe = glx+cnt+pointlist(1);
            yBottomStripe = gly+pointlist(2);
            
            [xSlide,ySlide]=polybool('union',xSlide,ySlide,xBottomStripe,yBottomStripe);
            
            cnt = cnt+0.2-(lw/nStripes);
        end%For
        
        
    case 'Pm'
        
        xMax =  5;
        xMin = -5;
        yMin = -0.5;
        
        xSlide = [xMin xMax xMax xMin];
        ySlide = [yMin-0.1 yMin-0.1 yMin-0.2 yMin-0.2];
        
        lw = 0.04; %linewidth of the lines below the ground
        glx = [0 -.1 lw-.1 lw]-0.5;
        gly = [-.5 -.6 -.6 -.5]-0.2;
        cnt = 0;
        %     for l = 1:11
        nStripes  = 6;
        for l = 1:nStripes
            
            
            xBottomStripe = glx+cnt+pointlist(1);
            yBottomStripe = gly+pointlist(2);
            
            [xSlide,ySlide]=polybool('union',xSlide,ySlide,xBottomStripe,yBottomStripe);
            
            cnt = cnt+0.2-(lw/nStripes);
        end%For
        
end

xOutline=polyUnionx;
yOutline=polyUniony;

% end

end%For




function [xOutline,yOutline,xSmallcircle,ySmallcircle] = getGroundOutline(DNA)
pointlist   = getGroundLocations(DNA);
spointlist  = size(pointlist);

xOutline = cell(spointlist(2));
yOutline = cell(spointlist(2));
xSmallcircle = cell(spointlist(2));
ySmallcircle = cell(spointlist(2));

for k = 1:spointlist(2)
    
    xTriangle = [-.5 .5 0]+pointlist(1,k);
    yTriangle = [-.5 -.5 0 ]+pointlist(2,k);
    
    t   = 0:0.1:2*pi;
    cx  = 0.1*cos(t) + pointlist(1,k);
    cy  = 0.1*sin(t) + pointlist(2,k);
    xSmallcircle{k} = 0.05*cos(t)+ pointlist(1,k);
    ySmallcircle{k} = 0.05*sin(t)+ pointlist(2,k);
            
    [polyUnionx,polyUniony]=polybool('union',xTriangle,yTriangle,cx,cy);
        
    lw = 0.04; %linewidth of the lines below the ground
    
    glx = [0 -.1 lw-.1 lw]-0.5;
    gly = [-.5 -.6 -.6 -.5];
    cnt = 0;
    for l = 1:11
        
        xBottomStripe = glx+cnt+pointlist(1,k);
        yBottomStripe = gly+pointlist(2,k);
        
        [polyUnionx,polyUniony]=polybool('union',polyUnionx,polyUniony,xBottomStripe,yBottomStripe);
        
        cnt = cnt+0.1-(lw/11);
    end%For
    
    xOutline{k}=polyUnionx;
    yOutline{k}=polyUniony;
    
end%For






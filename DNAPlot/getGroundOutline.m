



function [xOutline,yOutline,xSmallcircle,ySmallcircle,xlines,ylines] = getGroundOutline(DNA)
pointlist   = getGroundLocations(DNA);
spointlist  = size(pointlist);

xOutline = cell(1,spointlist(2));
yOutline = cell(1,spointlist(2));
xSmallcircle = cell(1,spointlist(2));
ySmallcircle = cell(1,spointlist(2));

% keyboard 

for k = 1:spointlist(2)
    
    w = 0.2; 
    h = 0.15;
    
    xTriangle = [-w/2  w/2 0]+pointlist(1,k);
    yTriangle = [-h -h 0 ]+pointlist(2,k);
    
    t   = 0:0.1:2*pi;
    cx  = 0.04*cos(t) + pointlist(1,k);
    cy  = 0.04*sin(t) + pointlist(2,k);
    xSmallcircle{k} = 0.04*cos(t)+ pointlist(1,k);
    ySmallcircle{k} = 0.04*sin(t)+ pointlist(2,k);
            
    [polyUnionx,polyUniony]=polybool('union',xTriangle,yTriangle,cx,cy);
    
    nLines = 5;
    hLines = h/5;
    wLines = hLines/tand(60);

    
    xlines{k} = [linspace(-w/2,w/2,nLines);linspace(-w/2,w/2,nLines)+wLines]+pointlist(1,k);
    ylines{k} = [-h*ones(1,nLines);(-h-hLines)*ones(1,nLines)]+pointlist(2,k);

    
    xOutline{k}=polyUnionx;
    yOutline{k}=polyUniony;
    
end%For






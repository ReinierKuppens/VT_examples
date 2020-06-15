%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %
%   Original Date: 29-06-2015                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p,pc,psl] = getGroundP(DNA,type)

p = cell(0);
pc = cell(0);
psl = cell(0);

[x,y,xs,ys,xSlide,ySlide] = groundCoords([0,0],type);

syms theta
R = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];
   
switch type 
    case 'P'
        for k = 1:numel(DNA.Ppar(1,:))
            
            XY = R(DNA.Ppar(3,k)+pi)*[x;y];
            XYslide    = R(DNA.Ppar(3,k)+pi)*[xSlide;ySlide];

            
            p{k}  = [XY(1,:)+DNA.Ppar(1,k); XY(2,:)+DNA.Ppar(2,k)];

            pc{k} = [xs+DNA.Ppar(1,k); ys+DNA.Ppar(2,k)];
            
            psl{k}  = [XYslide(1,:)+DNA.Ppar(1,k); XYslide(2,:)+DNA.Ppar(2,k)];

            
        end
    case 'Pm'
        for k = 1:numel(DNA.Pmpar(1,:))
            
            XY = R(DNA.Pmpar(3,k)+pi)*[x;y];
            p{k}  = [XY(1,:)+DNA.Pmpar(1,k); XY(2,:)+DNA.Pmpar(2,k)];
            pc{k} = [xs+DNA.Pmpar(1,k); ys+DNA.Pmpar(2,k)];  
% keyboard             
            XYslide    = R(DNA.Ppar(3,k)+pi)*[xSlide;ySlide];
            psl{k}  = [XYslide(1,:)+DNA.Ppar(1,k); XYslide(2,:)+DNA.Ppar(2,k)];

        end
        
end

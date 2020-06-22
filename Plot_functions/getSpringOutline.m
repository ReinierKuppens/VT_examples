%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %
%   Original Date: 29-06-2015                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x,y,nCoils,maxL,distUncoiled]=getSpringOutline(DNA,t,qss)

sSDNA = size(DNA.Spar);
nS = sSDNA(2);

if nS>0
    
    nCoils = zeros(1,nS);
    distUncoiled = zeros(1,nS);
    maxL = zeros(1,nS);
    x = cell(1,nS);
    y = cell(1,nS);
    
    for k = 1:nS
        
        p1          = DNA.Spar(1:2,k);
        p2          = DNA.Spar(3:4,k);
        Lspring     = getSpringLength(DNA,t,qss);
        nCoils(k)           = 15;
        distUncoiled(k)     = 0.1;
        maxL(k)             = max(Lspring(:,k))*2;
        nCoils(k)           = round(maxL(k)*4);

        [x{k}, y{k}] = springcoord(p1,p2,nCoils(k),maxL(k),distUncoiled(k));

    end
    
else
    x = cell(0);
    y = cell(0);
end




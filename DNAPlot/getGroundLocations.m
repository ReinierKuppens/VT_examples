%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %
%   Original Date: 29-06-2015                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [groundPointsFixed,groundPointsMove] = getGroundLocations(DNA)


incH = DNA2inc(DNA,'H');
incS = DNA2inc(DNA,'S');
incP = DNA2inc(DNA,'P');

groundPointsH = [];
groundPointsS = [];
groundPointsP = [];

if ~isempty(incH)
    groundPointsH = DNA.Hpar(1:2,incH(1,:)==1);
end
if ~isempty(incS)
    groundPointsS = DNA.Spar(1:2,incS(1,:)==1);
end
if ~isempty(incP)
    groundPointsP = DNA.Ppar(1:2,incP(1,:)==1);
end
% keyboard 
groundPointsFixed = [groundPointsH groundPointsS];
groundPointsMove  = [groundPointsP];

end
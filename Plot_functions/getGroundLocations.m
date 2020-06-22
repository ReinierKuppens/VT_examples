%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %
%   Original Date: 29-06-2015                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function groundPointsFixed = getGroundLocations(DNA)


incH = DNA2inc(DNA,'H');
incS = DNA2inc(DNA,'S');

groundPointsH = [];
groundPointsS = [];

if ~isempty(incH)
    groundPointsH = DNA.Hpar(1:2,incH(1,:)==1);
end
if ~isempty(incS)
    groundPointsS = DNA.Spar(1:2,incS(1,:)==1);
end

% keyboard 
groundPointsFixed = [groundPointsH groundPointsS];

end
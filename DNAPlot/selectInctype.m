%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %
%   Original Date: 29-06-2015                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function incstr = selectInctype(DNA,edgetype)

switch edgetype
    case 'H'
        incstr = DNA.incstr(DNA.edgelabel==1);
    case 'S'
        incstr = DNA.incstr(DNA.edgelabel==2);
    case 'P'
        incstr = DNA.incstr(DNA.edgelabel==3);
    case 'Pm'
        incstr = DNA.incstr(DNA.edgelabel==4);
    case 'all'
        incstr = DNA.incstr;
end%Switch
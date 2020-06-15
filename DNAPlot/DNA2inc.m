
% DNA2inc(DNA,edgetype) returns the incidence matrix that belongs to the
% specified edgetype It does so according to the pattern:
%
%     1     2     3     4     5     6     7     8     9    10   ...
%     -------------------------------------------------------
%     1     1     0     1     0     0     1     0     0     0   ...
%     1     0     1     0     1     0     0     1     0     0   ...
%     0     1     1     0     0     1     0     0     1     0   ...
%     0     0     0     1     1     1     0     0     0     1   ...
%     0     0     0     0     0     0     1     1     1     1   ...
%
%   Input:  DNA and edgetype, currently three options are possible:
%           'H'   -> for hinges
%           'S'   -> for springs
%           'R'   -> for gear ratios
%           'all' -> for complete DNA
%   Output: Incidence matrix
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %
%   Original Date: 30-06-2015                                             %
%   Last modified: 30-06-2015                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function inc = DNA2inc(DNA,edgetype)

edgelist    = selectEdgelist(DNA,edgetype);
sedgelist   = size(edgelist);

if ~isempty(edgelist)
    inc         = zeros(max(max(edgelist)),sedgelist(2));
    
    for k = 1:sedgelist(2)
        inc(edgelist(:,k),k) = 1;
    end%For
    
else
    inc = [];
end%If

end%Function


function edgelist = selectEdgelist(DNA,edgetype)

if      strcmp(edgetype,'H')
    edgelist = DNA2edgelist(DNA,'H');
    
elseif  strcmp(edgetype,'S')
    edgelist = DNA2edgelist(DNA,'S');
    
elseif strcmp(edgetype,'P')
    edgelist = DNA2edgelist(DNA,'P');
    
elseif strcmp(edgetype,'Pm')
    edgelist = DNA2edgelist(DNA,'Pm');
    
else
    edgelist = DNA2edgelist(DNA,'all');
end

end%Function













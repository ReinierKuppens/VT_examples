% DNA2edgelist(DNA,edgetype) returns the edgelist that is defined by the 
% incidence data of edgetype that are defined by the string. It does so 
% according to the pattern:
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
%           'all' -> for complete DNA
%   Output: 2xN vector containing the vertices that are connected in its
%           columns 
%
% % % % % DNA2edgelist(DNA,edgetype,data) returns the edgelist that is defined by 
% % % % % the incidence data of edgetype + aditional data of choice.
% % % % %
% % % % %   Input:  DNA,edgetype,data, currently three options are possible
% % % % %           'lab'       -> adds labels 
% % % % %           'par'       -> adds parameters    
% % % % %           'lab+par'   -> adds both labels and parameters
% % % % %
% % % % %   Output: depending on the input 'data' a differently sized matrix with
% % % % %   the data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %           
%   Original Date: 30-06-2015                                             %
%   Last modified: 30-06-2015                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edgelist = DNA2edgelist(DNA,edgetype)

incstr = selectInctype(DNA,edgetype);
if ~isempty(incstr)
    inc    = str2inc(incstr);
else
    edgelist = zeros(2,0);
    return
end

edgelist = zeros(2,length(incstr));

for k = 1:size(inc,2)
    edgelist(:,k) = find(inc(:,k)==1);
end

% if isempty(edgelist)
%     edgelist = zeros(0,2);
% end

end%Function





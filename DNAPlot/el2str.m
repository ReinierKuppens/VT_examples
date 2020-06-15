%str = el2str(el) gives the incidence string from an edgelist
%
%
%
%
%
% Author:        P.R. Kuppens 
% Date modified: 2017-12-20
%


function str = el2str(el)

el  = sort(el,1); 
str = 0.5.*(el(2,:) - 2).* (el(2,:)- 1) + el(1,:); 

%set str to zero for self loops 
str(el(1,:)==el(2,:))=0;

end
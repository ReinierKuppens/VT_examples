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

%Last version had the property that if el = [2;2] or [1;1]
% the incidence number is zero.
str(el(1,:)==el(2,:))=0;

%=====THIS IS THE OLD IMPLEMENTATIONS. IT IS VERY SLOW ====================
%     nE  = length(el(1,:));
%     nV  = max(max(el));
%     inc = zeros(nV,nE);
%     for  k = 1:nE
%         inc(el(:,k),k) = 1;
%     end
%     str = inc2str(inc);
%==========================================================================
end

% el = str2el(str) produces an edgelist from an incidence string. 
%
% The equations are based on the work done by Evert Kranendonk. 
%
%
%
% Author:        P.R. Kuppens 
% Date modified: 2018-01-24
%


function el = str2el(str)

el = zeros(2,length(str)); 

el(2,:) = floor(1.5 + 0.5.* sqrt( 1 + 8.*(str-1)));
el(1,:) = str - 0.5.*( el(2,:)-2) .* (el(2,:)-1 );

%=====THIS IS THE OLD IMPLEMENTATIONS. IT IS VERY SLOW ====================
% nE = length(str);
% el = zeros(2,nE);
% 
% inc = str2inc(str);
% 
% 
% for k = 1:size(inc,2)
%     el(:,k) = find(inc(:,k)==1);
% end
%==========================================================================
end
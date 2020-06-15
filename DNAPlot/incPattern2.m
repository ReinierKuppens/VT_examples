





%   M = GenCombinationMatrix(nCol) generates a matrix with all combinations
%   up to at least nCol columns. 
%   
%   [M,v] = GenCombinationMatrix(nCol) generates a matrix with all
%   combinations up to at least nE collumns and returns the minimum number
%   of vertices required.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %           
%   Original Date: 26-05-2015                                             %
%   Last modified: 26-05-2015                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = incPattern2(nCol)

nComb   = 0;        %Total number of columns in the combination matrix
c       = 1;        %Set counter to zero  

b = zeros(1,100);

while nComb<nCol    %While the Total number of columns is lower that the desired number
    c=c+1;          %Increase counter
    b(c) = nCr(c,2);
    nComb = nComb+(b(c)-b(c-1)); %Add collumns: (nCr(c,2)-nCr(c-1,2)) per vertices
end
v = c;          %Set counter equal to number of vertices

M = zeros(v,b(v));

for k = 1:(v-1)

    M(:,(b(k)+1:b(k+1))) = [eye(k);ones(1,k);zeros(v-k-1,k)];
    
end

 






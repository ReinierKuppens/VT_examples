% For inc = [1;0;0;0] it returns str = 0

function str = inc2str(mat)

if ~isempty(mat)
    
    sMat    = size(mat);
    v       = sMat(1);
    
    if sMat(1)==1 && sMat(2)==1 %for inc=1
        mat(2,1) = 0;
        v        = 2;
    end   
    
    M       = [];         %Empty combination Matrix
    c       = 1;          %Set counter to 1
    while c<v       %While the number of vertices processed is lower than the total number of vertices
        Mt = [];    %Clear Mt, which is a temporary matrix
        Mt = [eye(c);ones(1,c)];                    %Fill Mt with the structure
        s  = size(M);                               %Check size of M
        if isempty(M)                               %If M is an empty matrix
            M  = [M Mt] ;                           %It is possible to augment M and Mt directly
        else                                        %If M is not empty
            M  = [[M;zeros(1,s(2))] Mt];            %M and Mt are unequaly sized, so augment with appropriate number of zeros
        end
        c  = c+1;                                   %Increase the couter untill c=v
    end
    
    for k = 1:sMat(2)
        [~,str(k)]=ismember(mat(:,k).',M.','rows');
    end
    
else
    
    str = [];

end

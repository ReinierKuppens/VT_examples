



function Eg = getGravityEnergy(DNA,t,qss) 

g = 9.81;

for k = 1:numel(DNA.Mpar(:,1))-1
    for ii = 1:numel(t)
    
        m(k)     = DNA.Mpar(k+1,3);
        Eg(k,ii) = m(k)*g*qss(ii,k*3-2);
     
    end
end

Eg = Eg - min(Eg(:));






function plotPrisma(DNA)
hold on 

cdata       = [0.7 0.7 0.7];


[Pground,Pgroundc,Pslide]   = getGroundP(DNA,'P');
[Pmground,Pmgroundc,Pmslide] = getGroundP(DNA,'Pm');
 
nP  = numel(Pground);
nPm = numel(Pmground);
% keyboard 
for k = 1:nP
    Pgroundh{k}  = fill(Pground{k}(1,:),Pground{k}(2,:),cdata);
    Pgroundch{k} = fill(Pgroundc{k}(1,:),Pgroundc{k}(2,:),'w','Linewidth',2);
    Pslideh      = fill(Pslide{k}(1,:),Pslide{k}(2,:),cdata);
end
for k = 1:nPm
    Pmgroundh{k}  = fill(Pmground{k}(1,:),Pmground{k}(2,:),cdata);
    Pmgroundch{k} = fill(Pmgroundc{k}(1,:),Pmgroundc{k}(2,:),'w','Linewidth',2);
    Pmslideh      = fill(Pmslide{k}(1,:),Pmslide{k}(2,:),cdata);

end



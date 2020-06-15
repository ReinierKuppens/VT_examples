


function [] = plotMasses(DNA)


[massPoints,comPoint]   = getMassLocations(DNA);
nM                      = numel(massPoints);
% cdata                   = rand(nM,3);

if nM <= 3 
    cdata= [0.1       0.1       0.1;
            0.3467    0.5360    0.6907;
            0.9153    0.2816    0.2878];
else 
    cdata = linspecer(nM);
end
 
for k = 1:nM
    [x1{k},y1{k},x2{k},y2{k},x3{k},y3{k}] = getComOutline(comPoint{k},0.075); %ROW VECTORS
    [polx{k},poly{k},scx{k},scy{k}] = getMassOutline(massPoints{k}); %ROW VECTORS
end

for k =1:nM
    fill(polx{k},poly{k},cdata(k,:),'Edgecolor','k');
    alpha(0.5)
end

comColors = {'k','b','r','m','c'};
for k =1:nM
    for l =1:numel(scx{k})
        fill(scx{k}{l},scy{k}{l},'w','Linewidth',0.5);
    end
    
    fill(x3{k},y3{k},'w');
    plot(x3{k},y3{k},comColors{k});
    fill(x1{k},y1{k},comColors{k});
    fill(x2{k},y2{k},comColors{k});

end

end%Function










function cdata = getSpringColor(DNA,t,qss)

cdata       = cell(0,0);
[Energies]  = getEnergies(DNA,t,qss);
nS          = numel(Energies);
Es          = reshape([Energies{:}],length(t),nS);
maxE        = max(abs(Es(:)));

if ~isempty(Es)
    Es_normalized = Es./maxE;
end

nColor  = 256;
cmap    = jet(nColor);

Erank = abs(round(Es_normalized*nColor)); 
Erank(Erank==0)=1;

for ii = 1:length(t)
    for k = 1:nS
        cdata{k}(ii,:) = cmap(Erank(ii,k),:);
    end
end




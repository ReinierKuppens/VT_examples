


function hingeState = hingeState(DNA,t,qss)

[massPoints,comPoint]=getMassLocations(DNA);
nM = numel(massPoints);
R = @(angle) [cos(angle),-sin(angle);sin(angle),cos(angle)];

hingeState = [];
for k = 1:length(t)
    Hdata = [];
    
    for ii = 1:nM
        
        Rdata  = R(qss(k,ii*3))*[massPoints{ii}(:,1).'-comPoint{ii}(1);massPoints{ii}(:,2).'-comPoint{ii}(2)];
        RTdata = bsxfun(@plus,Rdata,[qss(k,ii*3-2);qss(k,ii*3-1)]);
        
        Hdata = [Hdata reshape(RTdata,1,numel(RTdata))];
        
    end
        hingeState = [hingeState;Hdata];
end
end







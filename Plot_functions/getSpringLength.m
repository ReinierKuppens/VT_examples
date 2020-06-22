

function Lspring = getSpringLength(DNA,t,qss)

  
sSDNA       = size(DNA.Spar);
nS          = sSDNA(2);
incS        = DNA2inc(DNA,'S');
comPoint    = getCompoints(DNA);
Lspring     = zeros(length(t),nS);

for k = 1:nS
    
    massnr  = find(incS(:,k)==1);    %Find entries in column
    ps1     = DNA.Spar(1:2,k);
    ps2     = DNA.Spar(3:4,k);
    
    for ii = 1:numel(t)
        
        if any(massnr==1)   %account for ground
            mnr = massnr(2)-1;
            mx  = qss(ii,mnr*3-2:mnr*3).';
            com = comPoint(mnr,:).';
            
            Lspring(ii,k) = norm(MGspringlength(ps1,ps2,[],[],com,mx));
        else
            mnr1 = massnr(1)-1;
            mnr2 = massnr(2)-1;
            
            mx   = [qss(ii,mnr1*3-2:mnr1*3) qss(ii,mnr2*3-2:mnr2*3)].';
            com1 = comPoint(mnr1,:).';
            com2 = comPoint(mnr2,:).';
            %             keyboard
            Lspring(ii,k) = norm(MMspringlength(ps1,ps2,[],[],com1,com2,mx));
        end
        
    end
    
end
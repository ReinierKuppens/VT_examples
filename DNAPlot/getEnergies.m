%[Es,dEs,Et,dEt] = getEnergies(DNA,varargin) give the total spring energy
%in the system. It gives it for both linear springs and torsion springs. In
%addition it gives the variation of the total spring energy. The total
%spring energy is defined as the sum of all individual springs under the
%condition of mobility. This means that if a spring is unable to move, or
%if a hinge is unable to rotate, it wont be included in the total energy.
%
% input: 1) DNA
%        2) h (optional) - sample time of numerical time integration
%        3) T (optional) - total duration of numerical integration
%
%   if 2) and 3) are not provided default values of h=0.05 and T=5 wil
%   used.
%
% output: Es  = total spring energy of linear springs over time
%         dEs = variation of total spring energy over time, if there's no
%         moving springs, variation of energy is inf
%         Et  = total torsional spring energy
%         dEt = variation of total torsion spring energy, if there's no
%         moving hinges energy is returned inf. Only for 0dof mechs. 


function Es = getEnergies(DNA,t,qss)

sSDNA       = size(DNA.Spar);
nS          = sSDNA(2);
incS        = DNA2inc(DNA,'S');
comPoint    = getCompoints(DNA);

Es = cell(1,nS);

for ii = 1:length(t)
    for k = 1:nS
        
        massnr = find(incS(:,k)==1);    %Find entries in column
        
        ps1 = DNA.Spar(1:2,k);
        ps2 = DNA.Spar(3:4,k);
        L0 = DNA.Spar(5,k);
        K  = DNA.Spar(6,k);
        
        if any(massnr==1)   %account for ground
           
            mnr = massnr(2)-1;
            com = comPoint(mnr,:).';
            
            MGEs = MGspringenergy(ps1,ps2,L0,K,com,[qss(ii,mnr*3-2);    %x1
                qss(ii,mnr*3-1);    %y1
                qss(ii,mnr*3-0)]);  %theta1
            
            Es{k}(ii) = MGEs;
            
        else
            
            mnr1 = massnr(1)-1;
            mnr2 = massnr(2)-1;
            com1 = comPoint(mnr1,:).';
            com2 = comPoint(mnr2,:).';
            
            
            MMEs = MMspringenergy(ps1,ps2,L0,K,com1,com2,[qss(ii,mnr1*3-2);     %x1
                qss(ii,mnr1*3-1);     %y1
                qss(ii,mnr1*3-0);     %theta1
                qss(ii,mnr2*3-2);     %x2
                qss(ii,mnr2*3-1);     %y2
                qss(ii,mnr2*3-0)]);   %theta2
            
            Es{k}(ii) = MMEs;
        end
    end
end
end










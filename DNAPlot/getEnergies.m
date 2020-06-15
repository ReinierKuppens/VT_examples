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


function [Es,dEs,Et,dEt] = getEnergies(DNA,varargin)

Es = [];
dEs = [];
Et = [];
dEt = [];

%Check if t and qss are provided
if isempty(varargin)
    h = 0.05;
    T = 5;
    [t,qss]=simulateDNA_mex(DNA,T,h);
else
    t   = varargin{1};
    qss = varargin{2};
end

sSDNA       = size(DNA.Spar);
sHDNA       = size(DNA.Hpar);
nS          = sSDNA(2);
nH          = sHDNA(2);
incS        = DNA2inc(DNA,'S');
incH        = DNA2inc(DNA,'H');
comPoint    = getCompoints(DNA);

% Es = zeros(length(t),1);
% Et = zeros(length(t),1);
%
% movingSprings = getMovingParts(incS,nS,qss); % THIS IS NOT RIGHT!!!!!!! This only checks if the two masses are moving or not. not the distance between the two points of the springs!
% movingHinges  = getMovingParts(incH,nH,qss);

% movingSprings = getMovingSprings(DNA,t,qss);
movingHinges  = getMovingHinges(DNA,t,qss);

movingSprings = 1:nS;


% keyboard 

Es = cell(1,numel(movingSprings));
Et = cell(1,numel(movingSprings));

% M1_prime         = cell(1,numel(alpha));
% M2_prime         = cell(1,numel(alpha));
% keyboard 
% 
for ii = 1:length(t)
    for k = movingSprings
        
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
            
%                         keyboard
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
% keyboard 

% dEs = gradient(Es);
for k = 1:movingSprings
    if sum(Es{k})>1e-6             %If there is energy at all, do
        dEs{k} = gradient(Es{k});
    else
        dEs{k} = inf(size(Es{k}));
    end
end

if nargout>2
    for ii = 1:length(t)
        for k = movingHinges
            
            massnr  = find(incH(:,k)==1);
            K       = DNA.Hpar(3,k);
            theta0  = DNA.Hpar(5,k);
            
            if any(massnr==1)
                mnr     = massnr(2)-1; %account for ground
                
                %Torsion Spring
                MGEt    = MGtorsionenergy(K,theta0,[qss(ii,mnr*3-2);    %x1
                    qss(ii,mnr*3-1);    %y1
                    qss(ii,mnr*3-0)]);  %theta
                Et(ii)  = Et(ii) + MGEt;
                
            else
                mnr1    = massnr(1)-1;
                mnr2    = massnr(2)-1;
                
                %Torsion spring
                MMEt    = MMtorsionenergy(K,theta0,[qss(ii,mnr1*3-2);   %x1
                    qss(ii,mnr1*3-1);   %y1
                    qss(ii,mnr1*3-0);   %theta1
                    qss(ii,mnr2*3-2);   %x2
                    qss(ii,mnr2*3-1);   %y2
                    qss(ii,mnr2*3-0)]); %theta2
                
                Et(ii)  = Et(ii) + MMEt;
                
            end
        end
    end
    %     keyboard
    %     dEt = gradient(Et);
    if sum(Et)>1e-6             %If there is energy at all, do
        dEt = gradient(Et);
    else
        dEt = inf(size(Et));
    end
    
end
end

%Check if springlength is changing, i.e. if spring move. 
function movingSprings = getMovingSprings(DNA,t,qss)

sSDNA       = size(DNA.Spar);
nS          = sSDNA(2);
incS        = DNA2inc(DNA,'S');
comPoint    = getCompoints(DNA);

movingSprings   = [];
Lspring         = zeros(length(t),1);

for k = 1:nS
    
    massnr  = find(incS(:,k)==1);    %Find entries in column
    ps1     = DNA.Spar(1:2,k);
    ps2     = DNA.Spar(3:4,k);
    
    for ii = 1:numel(t)
        
        if any(massnr==1)   %account for ground
            mnr = massnr(2)-1;
            
            mx  = qss(ii,mnr*3-2:mnr*3).';
            com = comPoint(mnr,:).';
            
            Lspring(ii) = norm(MGspringlength(ps1,ps2,[],[],com,mx));
        else
            mnr1 = massnr(1)-1;
            mnr2 = massnr(2)-1;
            
            mx   = [qss(ii,mnr1*3-2:mnr1*3) qss(ii,mnr2*3-2:mnr2*3)].';
            com1 = comPoint(mnr1,:).';
            com2 = comPoint(mnr2,:).';
            %             keyboard
            Lspring(ii) = norm(MMspringlength(ps1,ps2,[],[],com1,com2,mx));
        end
        
    end
    
    dLspring    = abs(diff(Lspring));
    sumdLspring = sum(dLspring);
    
%         keyboard
    %tolerance is dependent on stepsize h. For h = 0.05s nonmovable springs
    %have numerical variation of 4e-4m
    
    tol = 1;
    
    if sumdLspring>tol
        movingSprings = [movingSprings,k];
    end
    
    
end
end

%Check if the bodies rotate with respect to each other or ground.
function movingHinges = getMovingHinges(DNA,t,qss)

sHDNA       = size(DNA.Hpar);
nH          = sHDNA(2);
incH        = DNA2inc(DNA,'H');

movingHinges = [];

for k = 1:nH
    
    massnr = find(incH(:,k)==1);    %Find entries in column
    
    if any(massnr==1)   %account for ground
        mnr       = massnr(2)-1;
        theta1    = qss(:,mnr*3);
        dtheta1   = diff(theta1);
        sumdtheta = sum(abs(dtheta1));
        
%         keyboard 
        
    else
        mnr1 = massnr(1)-1;
        mnr2 = massnr(2)-1;
        
        theta1  = qss(:,mnr1*3);
        theta2  = qss(:,mnr2*3);
        
        theta12 = theta1-theta2;
        dtheta1 = diff(theta12);
        
        sumdtheta = sum(abs(dtheta1));
    end

%     sumdtheta
%     keyboard 
    tol = 1;

    if sumdtheta>tol
        movingHinges = [movingHinges,k];
    end
    
end
end













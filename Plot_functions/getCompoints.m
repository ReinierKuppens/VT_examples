%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %
%   Original Date: 29-06-2015                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function comPoint = getCompoints(DNA)

nM          = numel(DNA.Mpar(:,1))-1;    %Number of masses
comPoint    = zeros(nM,2);         %Initialize compoint array
inc         = DNA2inc(DNA,'all');

for k = 1:nM                                            %for all masses

    edgeIdx     = find(inc(k+1,:)==1);                     %verbonden edges per mass
    edgelabel   = DNA.edgelabel(edgeIdx);
     
    if numel(edgeIdx)==1                                %If there is only one edge
        points      = zeros(2);                         %Empty point array 
        points(1,:) = DNA.Mpar(k+1,1:2);                %Get first point out of mass data
        
        if edgelabel(1)==1            %for hinge               %If connection is hinge
            Hidx        = sum(DNA.edgelabel(1:edgeIdx(1))==edgelabel(1));
            points(2,:) = DNA.Hpar(1:2,Hidx).';
            
        elseif edgelabel(1)==2        %for spring
            Sidx = sum(DNA.edgelabel(1:edgeIdx(1))==edgelabel(1));
            vrtxIdx = find(inc(:,edgeIdx)==1);
            
            if k==vrtxIdx(1)-1
                points(2,:) = DNA.Spar(1:2,Sidx).';
                
            elseif k==vrtxIdx(2)-1
                points(2,:) = DNA.Spar(3:4,Sidx).';
                
            end
            
        elseif edgelabel(1)==3        %for prismatic joint
            Pidx = sum(DNA.edgelabel(1:edgeIdx(1))==edgelabel(1));
            points(2,:) = DNA.Ppar(1:2,Pidx).';
            
        elseif edgelabel(1)==4
            Pmidx = sum(DNA.edgelabel(1:edgeIdx(1))==edgelabel(1));
            points(2,:) = DNA.Pmpar(1:2,Pmidx).';
        end
        
        comPoint(k,:) = mean(points);
         
    else                                            %If there is more than one edge
        points = zeros(numel(edgeIdx),2);
        
        for ii = 1:numel(edgeIdx)
            if edgelabel(ii) == 1           %For hinge
                Hidx = sum(DNA.edgelabel(1:edgeIdx(ii))==edgelabel(ii));
                points(ii,:) = DNA.Hpar(1:2,Hidx).';
                
            elseif edgelabel(ii) == 2       %For spring
                vrtxIdx = find(inc(:,edgeIdx(ii))==1);
                Sidx = sum(DNA.edgelabel(1:edgeIdx(ii))==edgelabel(ii));

                if k == vrtxIdx(1)-1
                    points(ii,:) = DNA.Spar(1:2,Sidx).';
                    
                elseif k == vrtxIdx(2)-1
                    points(ii,:) = DNA.Spar(3:4,Sidx).';
                end
                
            elseif edgelabel(ii) == 3       %For prismatic joint
                 
                 Pidx = sum(DNA.edgelabel(1:edgeIdx(ii))==edgelabel(ii));
                 points(ii,:) = DNA.Ppar(1:2,Pidx).';
                 
            elseif edgelabel(ii) ==4 
                
                Pmidx = sum(DNA.edgelabel(1:edgeIdx(ii))==edgelabel(ii));
                points(ii,:) = DNA.Pmpar(1:2,Pmidx);
                
            end
        end
            
        comPoint(k,:) = mean(points);
        comPoint(k,:) = DNA.Mpar(k+1,1:2);
        
    end
end


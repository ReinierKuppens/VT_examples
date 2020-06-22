%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %
%   Original Date: 29-06-2015                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [points,comPoint] = getMassLocations(DNA)


nM          = numel(DNA.Mpar(:,1))-1;    %Number of masses
comPoint    = cell(1,nM);         %Initialize compoint array
inc         = DNA2inc(DNA,'all');

for k = 1:nM                                            %for all masses
    
    edgeIdx     = find(inc(k+1,:)==1);                     %verbonden edges per mass
    edgelabel   = DNA.edgelabel(edgeIdx);
    
    if isempty(edgeIdx)
        
        points{k}      = zeros(2);                         %Empty point array 
        points{k}(1,:) = DNA.Mpar(k+1,1:2);                %Get first point out of mass data    
        points{k}(2,:) = DNA.Mpar(k+1,1:2);                %Get first point out of mass data    
 
        comPoint{k} = mean(points{k});
        
    elseif numel(edgeIdx)==1                                %If there is only one edge
        points{k}      = zeros(2);                         %Empty point array 
        points{k}(1,:) = DNA.Mpar(k+1,1:2);                %Get first point out of mass data
        
        if edgelabel==1            %for hinge               %If connection is hinge
            
            Hidx            = sum(DNA.edgelabel(1:edgeIdx)==edgelabel);
            points{k}(2,:) = DNA.Hpar(1:2,Hidx).';
            
        elseif edgelabel==2        %for spring
            Sidx = sum(DNA.edgelabel(1:edgeIdx)==edgelabel);

            vrtxIdx = find(inc(:,edgeIdx)==1);
            if k==vrtxIdx(1)-1
                points{k}(2,:) = DNA.Spar(1:2,Sidx).';
            elseif k==vrtxIdx(2)-1
                points{k}(2,:) = DNA.Spar(3:4,Sidx).';
            end
        end

        comPoint{k} = mean(points{k});
 
        comPoint{k} = DNA.Mpar(k+1,1:2);

        
    else                                            %If there is more than one edge
        points{k} = zeros(numel(edgeIdx),2);
        
        for ii = 1:numel(edgeIdx)
            
            if edgelabel(ii) == 1           %For hinge
                Hidx = sum(DNA.edgelabel(1:edgeIdx(ii))==edgelabel(ii));
                points{k}(ii,:) = DNA.Hpar(1:2,Hidx).';
                
            elseif edgelabel(ii) == 2       %For spring
                vrtxIdx = find(inc(:,edgeIdx(ii))==1);
                Sidx = sum(DNA.edgelabel(1:edgeIdx(ii))==edgelabel(ii));

                if k == vrtxIdx(1)-1

                    points{k}(ii,:) = DNA.Spar(1:2,Sidx).';
                elseif k == vrtxIdx(2)-1

                    points{k}(ii,:) = DNA.Spar(3:4,Sidx).';
                end
                
            elseif edgelabel(ii) == 3       %For prismatic joint
                 Pidx = sum(DNA.edgelabel(1:edgeIdx(ii))==edgelabel(ii));
                 points{k}(ii,:) = DNA.Ppar(1:2,Pidx).';
            
            elseif edgelabel(ii) == 4
                 Pmidx = sum(DNA.edgelabel(1:edgeIdx(ii))==edgelabel(ii));
                 points{k}(ii,:) = DNA.Pmpar(1:2,Pmidx).';
                
            end
        end
        
%         keyboard 
        
        comPoint{k} = mean(points{k});
        
        comPoint{k} = DNA.Mpar(k+1,1:2);
         
        points{k} = [points{k};comPoint{k}];

%         points{k} = [points{k};];

        
    end
end
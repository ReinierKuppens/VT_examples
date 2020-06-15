

function [polyUnionx,polyUniony,smallcirclex,smallcircley] = getMassOutline(pointlist)
hold on

sortedPointlist = sortPolygonPoints(pointlist);         %Sort polygon points as to preven crossing lines
sizePointlist = size(pointlist);                        

radius1 = 0.1;          %Corner radius of each mass                  
radius2 = 0.04;         %Radius of each hinge axle
pointlist = [sortedPointlist; sortedPointlist(1,:)];

 
%Find all the tangent lines to the corner circles on each polygon vertex
polyUnionx = 0; 
polyUniony = 0; 
for k = 1:sizePointlist(1)      %For all corners points that for the polygon
    
    %Get normal vector
    p1   = pointlist(k,:);      
    p2   = pointlist(k+1,:);    
    diff = p1-p2;               
    dist = norm(diff);
    normalv    = diff./dist;
    
    A = [ 0  radius1;
         -radius1 0];
    
    %Translate each point such that it lies on the corner circle
    p3(:,k) = p1.'+A*normalv.';
    p4(:,k) = p1.'-A*normalv.';
    p5(:,k) = p2.'+A*normalv.';
    p6(:,k) = p2.'-A*normalv.';
    
    p = [p3(:,k) p4(:,k) p6(:,k) p5(:,k)];
    
    [polyUnionx,polyUniony] = poly2cw(polyUnionx,polyUniony);
    
    [p1all,p2all] = poly2cw(p(1,:),p(2,:));
    
    [polyUnionx, polyUniony] = polybool('union',polyUnionx,polyUniony,p1all,p2all);
     
end

    [polyUnionx,polyUniony] = poly2cw(polyUnionx,polyUniony);

    [polyUnionx,polyUniony] = polybool('union',polyUnionx,polyUniony,sortedPointlist(:,1),sortedPointlist(:,2));
        
%     polyUnionx(isnan(polyUnionx)==1)=[];
%     polyUniony(isnan(polyUniony)==1)=[];
%     
%     fill(polyUnionx,polyUniony,cdata,'Edgecolor','none')

idxcircle   = 0:0.1:2*pi;
for k = 1:sizePointlist(1)
    cx  = radius1*cos(idxcircle) + pointlist(k,1);
    cy  = radius1*sin(idxcircle) + pointlist(k,2);
    
    [cx,cy] = poly2cw(cx,cy);
    
    %This 'subtraction' caused invalid polygons in the next polybool with
    %'union' in some cases. Removing 'subtraction' fixes the error
%     [cxm, cym] = polybool('subtraction',cx,cy,polyUnionx,polyUniony);
    cxm =cx;cym=cy;
    
    cxm(isnan(cxm)==1)=[];
    cym(isnan(cym)==1)=[];

    [polyUnionx,polyUniony] = polybool('union',polyUnionx,polyUniony,cxm,cym);


end
alpha(0.3)

for k = 1:sizePointlist(1)
    smallcirclex{k} = radius2*cos(idxcircle)+ pointlist(k,1);
    smallcircley{k} = radius2*sin(idxcircle)+ pointlist(k,2);
end%For

    polyUnionx(isnan(polyUnionx)==1)=[];
    polyUniony(isnan(polyUniony)==1)=[];
    
    if size(polyUnionx,1)~=1
        polyUnionx = polyUnionx.';
        polyUniony = polyUniony.';
    end
   
k = convhull(polyUnionx,polyUniony);
polyUnionx = polyUnionx(k);
polyUniony = polyUniony(k); 
end%Function






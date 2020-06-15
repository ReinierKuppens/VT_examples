



function sortedPointlist = sortPolygonPoints(pointlist)

meanpoint = mean(pointlist,1); 

diff(:,1) = pointlist(:,1) - meanpoint(1);
diff(:,2) = pointlist(:,2) - meanpoint(2);

alpha = atan2(diff(:,2),diff(:,1));

[~,idx] = sort(alpha,'descend');

sortedPointlist = pointlist(idx,:);



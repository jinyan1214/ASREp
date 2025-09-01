function [xyz, WallElemNo, WallSegNo, SegElemNoList] = meshOnCurveWithFixNodes(cornersXYZ, wallCorners, ...
    appMeshSize, opening)

WallElemNo = zeros(size(wallCorners,1),1);
xyz = [];
SegElemNoList = [];
WallSegNo = zeros(size(WallElemNo));
for i = 1:size(wallCorners,1)
    wc1XYZ = cornersXYZ(wallCorners(i,1),:);
    wc2XYZ = cornersXYZ(wallCorners(i,2),:);
    endNode = 0;
    endNodeXYZ = wc1XYZ;
    if size(opening,1)==0
        openingOni = [];
    else
        openingOni1 = unique(opening(opening(:,1)==i,2));
        openingOni2 = unique(opening(opening(:,1)==i,4));
        openingOni = union(openingOni1,openingOni2);
    end
    
    
    for j = 1:length(openingOni)
            segElemNo = floor((openingOni(j)-endNode)/appMeshSize)+1;
            if floor((openingOni(j)-endNode)/appMeshSize) == ...
                    (openingOni(j)-endNode)/appMeshSize
               segElemNo =  floor((openingOni(j)-endNode)/appMeshSize);
            end
            WallElemNo(i) = WallElemNo(i) + segElemNo;
            headNodeXYZ = wc1XYZ + (wc2XYZ-wc1XYZ)/dist(wc1XYZ,wc2XYZ)*openingOni(j);
            for k=1:segElemNo
                xyz = [xyz; endNodeXYZ+(headNodeXYZ-endNodeXYZ)/segElemNo*k];
            end
            endNodeXYZ = xyz(end,:);
            endNode = openingOni(j);
            WallSegNo(i) = WallSegNo(i)+1;
            SegElemNoList = [SegElemNoList; segElemNo];
    end
    segElemNo = floor((norm(wc2XYZ-wc1XYZ)-endNode)/appMeshSize)+1;
    if floor((norm(wc2XYZ-wc1XYZ)-endNode)/appMeshSize) == ...
            (norm(wc2XYZ-wc1XYZ)-endNode)/appMeshSize
       segElemNo =  floor((norm(wc2XYZ-wc1XYZ)-endNode)/appMeshSize);
    end
    WallElemNo(i) = WallElemNo(i) + segElemNo;
    headNodeXYZ = wc2XYZ;
    for k=1:segElemNo
        xyz = [xyz; endNodeXYZ+(headNodeXYZ-endNodeXYZ)/segElemNo*k];
    end
    WallSegNo(i) = WallSegNo(i)+1;
    SegElemNoList = [SegElemNoList; segElemNo];
end

end
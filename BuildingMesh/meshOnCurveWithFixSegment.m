function xyz = meshOnCurveWithFixSegment(cornersXYZ, exterCornerXYZ, wallCorners, ...
    SegElemNoList, opening)
xyz = [];
wallSegAcc = 0;
for i = 1:size(wallCorners,1)
    wc1XYZ = cornersXYZ(wallCorners(i,1),:);
    wc2XYZ = cornersXYZ(wallCorners(i,2),:);
    endNode = 0;
    endNodeXYZ = wc1XYZ;
    
    if size(opening,1) == 0
        openingOniLocal = [];
        openingOni=[];
    else
        openingOni1 = unique(opening(opening(:,1)==i,2));
        openingOni2 = unique(opening(opening(:,1)==i,4));
        openingOni = union(openingOni1,openingOni2);
        openingOniLocal = zeros(size(openingOni));
    end
    
    exter1XYZ = exterCornerXYZ(wallCorners(i,1),:);
    exter2XYZ = exterCornerXYZ(wallCorners(i,2),:);
    for j = 1:length(openingOni)
        aVect = wc1XYZ-exter1XYZ;
        bVect = (wc2XYZ-wc1XYZ)/norm(wc2XYZ-wc1XYZ);
        projLength = (aVect*(bVect.'));
        openingOniLocal(j) = openingOni(j)-projLength;
    end
    walliSegNo = length(openingOni) + 1;
    for j = 1:length(openingOni)
            segElemNo = SegElemNoList(wallSegAcc+j);           
            headNodeXYZ = wc1XYZ + (wc2XYZ-wc1XYZ)/dist(wc1XYZ,wc2XYZ)*openingOniLocal(j);
            for k=1:segElemNo
                xyz = [xyz; endNodeXYZ+(headNodeXYZ-endNodeXYZ)/segElemNo*k];
            end
            endNodeXYZ = xyz(end,:);
            endNode = openingOniLocal(j);
    end
    segElemNo = SegElemNoList(wallSegAcc + length(openingOni)+1);
    headNodeXYZ = wc2XYZ;
    for k=1:segElemNo
        xyz = [xyz; endNodeXYZ+(headNodeXYZ-endNodeXYZ)/segElemNo*k];
    end
    wallSegAcc = wallSegAcc + walliSegNo;
end

end
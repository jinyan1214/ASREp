function xyz = meshOnCurve(cornersXYZ, wallCorners, noOfElems)
xyz = zeros(sum(noOfElems),3);
nodeDefined = 0;
for i = 1:length(noOfElems)
   corner1XYZ = cornersXYZ(wallCorners(i,1),:);
   corner2XYZ = cornersXYZ(wallCorners(i,2),:);
   for j = 1:noOfElems(i)
       xyz(nodeDefined + j,:) = corner1XYZ + ...
           (corner2XYZ-corner1XYZ)/noOfElems(i)*j;
   end
   nodeDefined = nodeDefined + noOfElems(i);
end

end
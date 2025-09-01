%(cornersXYZ, wallCorners, offset)
function [newCornersXYZ, newNodes] = offsetCurveOut(cornersXYZ, wallCorners, offset, wallnodes, wallNoElem)
    newCornersXYZ = zeros(size(cornersXYZ));
    newNodes = zeros(size(wallnodes));
    for i = 1:size(wallCorners,1)
       w1c1 = cornersXYZ(wallCorners(i,1),:);
       w1c2 = cornersXYZ(wallCorners(i,2),:);
       offset1 = offset(i);
       if i<size(wallCorners,1)
            w2c1 = cornersXYZ(wallCorners(i+1,1),:);
            w2c2 = cornersXYZ(wallCorners(i+1,2),:);
            offset2 = offset(i+1);
       else
            w2c1 = cornersXYZ(wallCorners(1,1),:);
            w2c2 = cornersXYZ(wallCorners(1,2),:);
            offset2 = offset(1);
       end
%        theta = pi/2;
%        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
       R = [0 1 0; -1 0 0;0 0 1];%clockwise 90 degree
       w1 = w1c2 - w1c1;
       w1Perpen = R*(w1.');
       w1c1Off = w1c1 + w1Perpen'/norm(w1Perpen)*offset1;
       w1c2Off = w1c2 + w1Perpen'/norm(w1Perpen)*offset1;
       
       w2 = w2c2 - w2c1;
       w2Perpen = R*(w2.');
       w2c1Off = w2c1 + w2Perpen'/norm(w2Perpen)*offset2;
       w2c2Off = w2c2 + w2Perpen'/norm(w2Perpen)*offset2;
       
       newCorner = intersect2lineGiven4Points2D(...
           w1c1Off, w1c2Off, w2c1Off, w2c2Off);
       
       exterNo = sum(wallNoElem(1:i-1))+1;
       
       if i<size(wallCorners,1)
           newCornersXYZ(i+1,:) = [newCorner, w1c1(3)];
           
       else
           newCornersXYZ(1,:) = [newCorner, w1c1(3)];
       end
       
       newNodes(exterNo:exterNo+wallNoElem(i),:) = wallnodes(exterNo:exterNo+wallNoElem(i)) + ...
                w1Perpen.'/norm(w1Perpen)*offset(i);
       newNodes(exterNo+wallNoElem(i),:) = [newCorner, w1c1(3)];
       
       
    end

end
%(cornersXYZ, wallCorners, offset)
function newCornersXYZ = offsetCurveIn(cornersXYZ, wallCorners, offset)
    newCornersXYZ = zeros(size(cornersXYZ));
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
       R = [0 -1; 1 0];
       w1 = w1c2 - w1c1;
       w1Perpen = R*(w1(1:2).');
       w1c1Off = w1c1(1:2) + w1Perpen'/norm(w1Perpen)*offset1;
       w1c2Off = w1c2(1:2) + w1Perpen'/norm(w1Perpen)*offset1;
       if size(wallCorners,1)==1
           newCornersXYZ(1,1:2) = w1c1Off;
           newCornersXYZ(2,1:2) = w1c2Off;
           return
       end
       
       w2 = w2c2 - w2c1;
       w2Perpen = R*(w2(1:2).');
       w2c1Off = w2c1(1:2) + w2Perpen'/norm(w2Perpen)*offset2;
       w2c2Off = w2c2(1:2) + w2Perpen'/norm(w2Perpen)*offset2;
       
       [xi,yi] = linexline([w1c1Off(1), w1c2Off(1)], [w1c1Off(2), w1c2Off(2)],...
           [w2c1Off(1), w2c2Off(1)], [w2c1Off(2), w2c2Off(2)], 0);
       newCorner = [xi, yi];
       if i<size(wallCorners,1)
           newCornersXYZ(i+1,1:2) = newCorner;
       else
           newCornersXYZ(1,1:2) = newCorner;
       end
    end

end
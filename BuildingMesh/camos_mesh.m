% demo3 of MasonMesh tool (Add partition walls to facade)
% This demo create a mesh for a 48x16m masonary facade sitting on a 1m
% thick foundation. 
% The facade walls have openings as defined in "Camós, C., Molins, C., 
% & Arnau, O. (2014). Case study of damage on masonry buildings
% produced by tunneling induced settlements. International Journal of 
% Architectural Heritage, 8(4), 602-625."
%
% Revision history:
%   Jinyan Zhao, jinyan_zhao@berkeley.edu, Jan 2023
% Cite As
%   
% 
clear 
close all
%% Input 
% Define the coordinate of facade corners
baseCornersXYZ = [ -24, -8, 0;
                  24, -8, 0;
                  24,  8, 0;
                 -24,  8, 0];
% Define the connectivity of the corners, the define of base wall has 
% to be counter clock wise             
baseWall = [1 2;
            2 3;
            3 4;
            4 1];
% Thickness of each wall, each wall can have different thickness       
WallThick = 0.2*ones(4,1);
% Thickness of foundation below each wall, all foundations can have
% different width
foundationWidth = 0.3*ones(size(baseWall,1),1);
% Height of the walls
WallHeight = 3;
% Height of the foundation
foundationHeight = 0.25;
% Thickness of partion wall
partWallThick = 0.04;
%--------------------Meshing parameters------------------------------------
WallThickNoElem = 2; % Number of elements perpendicular to the facade plain
WallExterElemApproSize = 1;% Approximate element size in the facade plain
heightElemSize = 1; % Approximage facade element size in the vertical direction
foundationOutNoElem = 1; % Number of foundation element outside of the facade 
foundationInNoElem = 1; % Number of foundation element inside of the facade
foundationHeighNoElem = 2; % Number of foundation element in the vertical direction
partWallThickNoElem = 1; % Number of elements perpendicular to the partition plain
%--------------------- input of opening------------------------------------
% opening is defined as wall_no, x1/y1 z1 x2/y2 z2 
opening = [];
opening = [opening; [1 1.3 1 2.3 2]];
opening = [opening; [1 3.6 0 4.4 2]];
opening = [opening; [1 5.7 1 6.7 2]];
opening = [opening; [1 9.3 1 10.3 2]];
opening = [opening; [1 11.6 0 12.4 2]];
opening = [opening; [1 13.7 1 14.7 2]];
opening = [opening; [1 17.3 1 18.3 2]];
opening = [opening; [1 19.6 0 20.4 2]];
opening = [opening; [1 21.7 1 22.7 2]];
opening2 = opening;
opening2(:,2) = opening2(:,2)+24;
opening2(:,4) = opening2(:,4)+24;
opening = [opening; opening2];
opening2 = opening;
opening2(:,1) = 3;
opening = [opening; opening2];
opening = [opening; [1 8-partWallThick/2 0 8+partWallThick/2 0]];
opening = [opening; [1 16-partWallThick/2 0 16+partWallThick/2 0]];
opening = [opening; [1 24-partWallThick/2 0 24+partWallThick/2 0]];
opening = [opening; [1 32-partWallThick/2 0 32+partWallThick/2 0]];
opening = [opening; [1 40-partWallThick/2 0 40+partWallThick/2 0]];
opening = [opening; [3 8-partWallThick/2 0 8+partWallThick/2 0]];
opening = [opening; [3 16-partWallThick/2 0 16+partWallThick/2 0]];
opening = [opening; [3 24-partWallThick/2 0 24+partWallThick/2 0]];
opening = [opening; [3 32-partWallThick/2 0 32+partWallThick/2 0]];
opening = [opening; [3 40-partWallThick/2 0 40+partWallThick/2 0]];
opening = [opening; [2 8-partWallThick/2 0 8+partWallThick/2 0]];
opening = [opening; [4 8-partWallThick/2 0 8+partWallThick/2 0]];
% split window and door into 2 elems
addNodes = [1.8, 4, 6.2, 9.8, 12, 14.2, 17.8, 20, 22.2, 25.8, 28, 30.2, 33.8, 36, 38.2, 41.8, 44, 46.2];
for i  = 1:length(addNodes)
    opening = [opening; [1 addNodes(i) 0 addNodes(i) 0]];
end
for i  = 1:length(addNodes)
    opening = [opening; [3 addNodes(i) 0 addNodes(i) 0]];
end
%% Generate facade nodes (wallNodesXYZ) and mesh (elem2n)
% nodes
WallLength = zeros(size(WallThick));
WallOpenNo = zeros(size(WallThick));
[baseExterNodesXYZ, WallNoElem, wallSegNo, segElemNo]= meshOnCurveWithFixNodes(baseCornersXYZ, baseWall, ...
    WallExterElemApproSize, opening);
baseNodesXYZ = zeros(sum(WallNoElem)*(WallThickNoElem+1),3);
baseNodesXYZ(1:sum(WallNoElem),:) = baseExterNodesXYZ;
outerCornersXYZ = baseCornersXYZ;
for i=1:WallThickNoElem
    % find corners
    newCornersXYZ = offsetCurveIn(outerCornersXYZ, baseWall, WallThick/WallThickNoElem);
    % find mesh on curve
    baseNodesXYZ(1+sum(WallNoElem)*i:sum(WallNoElem)*(i+1),:) = ...
        meshOnCurveWithFixSegment(newCornersXYZ, baseCornersXYZ, baseWall, ...
            segElemNo, opening);   
    outerCornersXYZ = newCornersXYZ;
end
baseInnerCornersXYZ = newCornersXYZ;
if size(opening,1)>0
    openingZ = union(unique(opening(:,3)),unique(opening(:,5)));
    openingZ = [0; openingZ; WallHeight];
    openingZ = unique(openingZ);
else
    openingZ = [0; WallHeight];
end
endHeight = 0;
heightZ = 0;
for i=2:length(openingZ)
   headHeight = openingZ(i);
   segElemNo = floor((headHeight-endHeight)/heightElemSize)+1;
   if floor((headHeight-endHeight)/heightElemSize) == (headHeight-endHeight)/heightElemSize
       segElemNo = floor((headHeight-endHeight)/heightElemSize);
   end
   
   for k = 1:segElemNo
      heightZ = [heightZ; endHeight + k*(headHeight-endHeight)/segElemNo]; 
   end
   endHeight = heightZ(end);
end
% Elements
elemPerLayerWall = sum(WallNoElem)*WallThickNoElem;
WallHeightNoElem = length(heightZ)-1;
elem2n = zeros(elemPerLayerWall*WallHeightNoElem, 8);
for j = 1:WallThickNoElem
    for i=1:sum(WallNoElem)
       elem2n(i+(j-1)*sum(WallNoElem),1) = i + (j-1)*sum(WallNoElem);
       elem2n(i+(j-1)*sum(WallNoElem),2) = i+1 + (j-1)*sum(WallNoElem);
       elem2n(i+(j-1)*sum(WallNoElem),3) = i+1 + sum(WallNoElem) + (j-1)*sum(WallNoElem);
       elem2n(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem);
       if i==sum(WallNoElem)
           elem2n(i+(j-1)*sum(WallNoElem),1) = i + (j-1)*sum(WallNoElem);
           elem2n(i+(j-1)*sum(WallNoElem),2) = 1 + (j-1)*sum(WallNoElem);
           elem2n(i+(j-1)*sum(WallNoElem),3) = 1 + sum(WallNoElem) + (j-1)*sum(WallNoElem);
           elem2n(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem);
       end
    end
end
elem2n(1:elemPerLayerWall,5:8) = ...
    elem2n(1:elemPerLayerWall,1:4) + size(baseNodesXYZ,1);
wallNodesXYZ = zeros(size(baseNodesXYZ,1)*(WallHeightNoElem+1),3);
wallNodesXYZ(1:size(baseNodesXYZ,1),:) = baseNodesXYZ;
for i = 1:WallHeightNoElem-1
   elem2n(i*elemPerLayerWall+1:(i+1)*elemPerLayerWall,:) = ...
       elem2n((i-1)*elemPerLayerWall+1:i*elemPerLayerWall,:) + ...
       size(baseNodesXYZ,1);
   wallNodesXYZ(i*size(baseNodesXYZ,1)+1:(i+1)*size(baseNodesXYZ,1),1:2) =...
       baseNodesXYZ(:,1:2);
   wallNodesXYZ(i*size(baseNodesXYZ,1)+1:(i+1)*size(baseNodesXYZ,1),3) = ...
       heightZ(i+1);
end
wallNodesXYZ(WallHeightNoElem*size(baseNodesXYZ,1)+1:(WallHeightNoElem+1)*size(baseNodesXYZ,1),1:2) =...
       baseNodesXYZ(:,1:2);
wallNodesXYZ(WallHeightNoElem*size(baseNodesXYZ,1)+1:(WallHeightNoElem+1)*size(baseNodesXYZ,1),3) = ...
       heightZ(end);
PlotMesh(wallNodesXYZ,elem2n, 0)
%% Remove open element
windowElemNo = [];
for i = 1:size(opening,1)
    wall_no = opening(i, 1);
    if wall_no>1
       exterNoX = sum(WallNoElem(1:(wall_no-1)))+1;
       innerNoX = sum(WallNoElem(1:(wall_no-1)))+1;
    else
       exterNoX = 0+1;
       innerNoX = 0+1;
    end
    wc1xyz = baseCornersXYZ(baseWall(wall_no,1),:);
    wc2xyz = baseCornersXYZ(baseWall(wall_no,2),:);
    for j = exterNoX:sum(WallNoElem(1:wall_no))-1
        for k = 1:WallHeightNoElem
            for l = 1:WallThickNoElem
                elem_no = j + (k-1)*sum(WallNoElem)*WallThickNoElem + ...
                    (l-1)*sum(WallNoElem);
                meanX = mean(wallNodesXYZ(elem2n(elem_no,:),1));
                meanY = mean(wallNodesXYZ(elem2n(elem_no,:),2));
                meanZ = mean(wallNodesXYZ(elem2n(elem_no,:),3));
                projXYZ = projectOnLine([meanX,meanY,meanZ], wc1xyz, wc2xyz);
                projXYZ(3)=0;
                if opening(i,2)<dist(projXYZ, wc1xyz) && ...
                        opening(i,4)>dist(projXYZ, wc1xyz) && ...
                        opening(i,3)<meanZ && ...
                        opening(i,5)>meanZ
                    windowElemNo = [windowElemNo,elem_no];
                end
            end 
        end
    end       
end
elem2n(windowElemNo,:) = [];
PlotMesh(wallNodesXYZ,elem2n, 0)
%% Add foundation nodes and element to get wholeNodesXYZ and wholeElem2n
startFoundationElem = length(elem2n)+1;
% nodes
foundationNodesXYZ = zeros(sum(WallNoElem)*(WallThickNoElem+1+foundationOutNoElem+foundationInNoElem) * ...
    (foundationHeighNoElem+1) - size(baseNodesXYZ,1),3);
% foundationNodesXYZ(1:sum(WallNoElem),:) = offsetWallNodes(wallNodesXYZ(1:sum(WallNoElem),:),...
%     wallNodesXYZ(1+sum(WallNoElem):sum(WallNoElem)*2,:),WallThick/WallThickNoElem,...
%     (foundationWidth-WallThick)/2, WallNoElem);
foundationNodesXYZ(1:sum(WallNoElem),:) = offsetWallNodes(wallNodesXYZ(1:sum(WallNoElem),:),...
    wallNodesXYZ(1+sum(WallNoElem):sum(WallNoElem)*2,:),WallThick/WallThickNoElem,...
    (foundationWidth-WallThick)/2);

for i = 1:foundationOutNoElem-1
%     foundationNodesXYZ(1+i*sum(WallNoElem):(i+1)*sum(WallNoElem),:) = ...
%         offsetWallNodes(wallNodesXYZ(1:sum(WallNoElem),:),...
%         wallNodesXYZ(1+sum(WallNoElem):sum(WallNoElem)*2,:),WallThick/WallThickNoElem,...
%         (foundationWidth-WallThick)/2-i*(foundationWidth-WallThick)/2/foundationOutNoElem, WallNoElem);
    foundationNodesXYZ(1+i*sum(WallNoElem):(i+1)*sum(WallNoElem),:) = ...
        offsetWallNodes(wallNodesXYZ(1:sum(WallNoElem),:),...
        wallNodesXYZ(1+sum(WallNoElem):sum(WallNoElem)*2,:),WallThick/WallThickNoElem,...
        (foundationWidth-WallThick)/2-i*(foundationWidth-WallThick)/2/foundationOutNoElem);
end
for i = 1:foundationInNoElem
%     foundationNodesXYZ(1+(i-1+foundationOutNoElem)*sum(WallNoElem):(i-1+1+foundationOutNoElem)*sum(WallNoElem),:) = ...
%         offsetWallNodes(wallNodesXYZ(1+sum(WallNoElem)*WallThickNoElem:sum(WallNoElem)+sum(WallNoElem)*WallThickNoElem,:),...
%         wallNodesXYZ(1+sum(WallNoElem)*(WallThickNoElem-1):sum(WallNoElem)*WallThickNoElem,:),WallThick/WallThickNoElem,...
%         (foundationWidth-WallThick)/2-(i-1)*(foundationWidth-WallThick)/2/foundationInNoElem, WallNoElem);
    foundationNodesXYZ(1+(i-1+foundationOutNoElem)*sum(WallNoElem):(i-1+1+foundationOutNoElem)*sum(WallNoElem),:) = ...
        offsetWallNodes(wallNodesXYZ(1+sum(WallNoElem)*WallThickNoElem:sum(WallNoElem)+sum(WallNoElem)*WallThickNoElem,:),...
        wallNodesXYZ(1+sum(WallNoElem)*(WallThickNoElem-1):sum(WallNoElem)*WallThickNoElem,:),WallThick/WallThickNoElem,...
        (foundationWidth-WallThick)/2-(i-1)*(foundationWidth-WallThick)/2/foundationInNoElem);
end

foundationNodesXYZ(sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)+1:...
    sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)+sum(WallNoElem)*foundationOutNoElem,:) = ...
    foundationNodesXYZ(1:sum(WallNoElem)*(0+foundationOutNoElem),:);

foundationNodesXYZ(sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
    +sum(WallNoElem)*foundationOutNoElem+1:...
    sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
    +sum(WallNoElem)*foundationOutNoElem + size(baseNodesXYZ,1),:) = ...
    baseNodesXYZ;

foundationNodesXYZ(sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
    +sum(WallNoElem)*foundationOutNoElem+size(baseNodesXYZ,1)+1:...
    sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
    +sum(WallNoElem)*foundationOutNoElem + size(baseNodesXYZ,1)...
    +sum(WallNoElem)*foundationInNoElem,:) = ...
    foundationNodesXYZ(sum(WallNoElem)*foundationOutNoElem+1:...
    sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem),:);


foundationNodesXYZ(sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)+1:...
    sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)*2 + size(baseNodesXYZ,1),3) = ...
    0 - foundationHeight/foundationHeighNoElem;

foundationLayerNoNodes = sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem) + size(baseNodesXYZ,1);
for i = 1:foundationHeighNoElem-1
    foundationNodesXYZ(sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
        +i*foundationLayerNoNodes + 1:...
        sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
        +(i+1)*foundationLayerNoNodes,:) = ...
    foundationNodesXYZ(sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
        +(i-1)*foundationLayerNoNodes + 1:...
        sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
        +(i)*foundationLayerNoNodes,:);

    foundationNodesXYZ(sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
        +i*foundationLayerNoNodes + 1:...
        sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
        +(i+1)*foundationLayerNoNodes,3) = 0 - foundationHeight/foundationHeighNoElem*(i+1);
end

wholeNodesXYZ = [wallNodesXYZ;foundationNodesXYZ];

elem2nFound = zeros((WallThickNoElem+foundationOutNoElem+foundationInNoElem)*...
    sum(WallNoElem) * foundationHeighNoElem,8);
for j= 1:foundationOutNoElem
    for i = 1:sum(WallNoElem)
       elem2nFound(i+(j-1)*sum(WallNoElem),1) = i+(j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
       elem2nFound(i+(j-1)*sum(WallNoElem),2) = i+1 + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
       elem2nFound(i+(j-1)*sum(WallNoElem),3) = i+1 + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
       elem2nFound(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
       elem2nFound(i+(j-1)*sum(WallNoElem),5) = i+(j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),6) = i+1 + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),7) = i+1 + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       if j==foundationOutNoElem
           elem2nFound(i+(j-1)*sum(WallNoElem),3) = i+1;
           elem2nFound(i+(j-1)*sum(WallNoElem),4) = i;
           elem2nFound(i+(j-1)*sum(WallNoElem),7) = i+1 + size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
               + sum(WallNoElem)*foundationOutNoElem;
           elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
               + sum(WallNoElem)*foundationOutNoElem;
       end
       if i==sum(WallNoElem)
           elem2nFound(i+(j-1)*sum(WallNoElem),1) = i + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
           elem2nFound(i+(j-1)*sum(WallNoElem),2) = 1 + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
           elem2nFound(i+(j-1)*sum(WallNoElem),3) = 1 + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
           elem2nFound(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
           elem2nFound(i+(j-1)*sum(WallNoElem),5) = i + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),6) = 1 + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),7) = 1 + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           if j==foundationOutNoElem
               elem2nFound(i+(j-1)*sum(WallNoElem),3) = 1;
               elem2nFound(i+(j-1)*sum(WallNoElem),4) = i;
               elem2nFound(i+(j-1)*sum(WallNoElem),7) = 1 + size(wallNodesXYZ,1)...
                   + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
                   + sum(WallNoElem)*foundationOutNoElem;
               elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + size(wallNodesXYZ,1)...
                   + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
                   + sum(WallNoElem)*foundationOutNoElem;
           end
       end  
    end
end

for j= 1+foundationOutNoElem:WallThickNoElem+foundationOutNoElem
    for i = 1:sum(WallNoElem)
       elem2nFound(i+(j-1)*sum(WallNoElem),1) = i+(j-foundationOutNoElem-1)*sum(WallNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),2) = i+1 + (j-foundationOutNoElem-1)*sum(WallNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),3) = i+1 + sum(WallNoElem) + (j-foundationOutNoElem-1)*sum(WallNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-foundationOutNoElem-1)*sum(WallNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),5) = i+(j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),6) = i+1 + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),7) = i+1 + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       if i==sum(WallNoElem)
           elem2nFound(i+(j-1)*sum(WallNoElem),1) = i + (j-foundationOutNoElem-1)*sum(WallNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),2) = 1 + (j-foundationOutNoElem-1)*sum(WallNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),3) = 1 + sum(WallNoElem) + (j-foundationOutNoElem-1)*sum(WallNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-foundationOutNoElem-1)*sum(WallNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),5) = i + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),6) = 1 + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),7) = 1 + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       end  
    end
end

for j= 1+foundationOutNoElem+WallThickNoElem:...
        WallThickNoElem+foundationOutNoElem+foundationInNoElem
    for i = 1:sum(WallNoElem)
       elem2nFound(i+(j-1)*sum(WallNoElem),1) = i+(j-WallThickNoElem-1-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
       elem2nFound(i+(j-1)*sum(WallNoElem),2) = i+1 + (j-WallThickNoElem-1-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
       elem2nFound(i+(j-1)*sum(WallNoElem),3) = i+1 + sum(WallNoElem) + (j-WallThickNoElem-1-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
       elem2nFound(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-WallThickNoElem-1-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
       elem2nFound(i+(j-1)*sum(WallNoElem),5) = i+(j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),6) = i+1 + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),7) = i+1 + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       if j==1+foundationOutNoElem+WallThickNoElem
           elem2nFound(i+(j-1)*sum(WallNoElem),1) = i + sum(WallNoElem)*WallThickNoElem;
           elem2nFound(i+(j-1)*sum(WallNoElem),2) = i + 1 + sum(WallNoElem)*WallThickNoElem;
           elem2nFound(i+(j-1)*sum(WallNoElem),5) = i + size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
               + sum(WallNoElem)*(foundationOutNoElem + WallThickNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),6) = i + 1 + size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
               + sum(WallNoElem)*(foundationOutNoElem + WallThickNoElem);
       end
       if i==sum(WallNoElem)
           elem2nFound(i+(j-1)*sum(WallNoElem),1) = i+(j-WallThickNoElem-1-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
           elem2nFound(i+(j-1)*sum(WallNoElem),2) = 1 + (j-WallThickNoElem-1-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
           elem2nFound(i+(j-1)*sum(WallNoElem),3) = 1 + sum(WallNoElem) + (j-WallThickNoElem-1-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
           elem2nFound(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-WallThickNoElem-1-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
           elem2nFound(i+(j-1)*sum(WallNoElem),5) = i+(j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),6) = 1 + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),7) = 1 + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           if j==1+foundationOutNoElem+WallThickNoElem
               elem2nFound(i+(j-1)*sum(WallNoElem),1) = i + sum(WallNoElem)*WallThickNoElem;
               elem2nFound(i+(j-1)*sum(WallNoElem),2) = 1 + sum(WallNoElem)*WallThickNoElem;
               elem2nFound(i+(j-1)*sum(WallNoElem),5) = i + size(wallNodesXYZ,1)...
                   + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
                   + sum(WallNoElem)*(foundationOutNoElem + WallThickNoElem);
               elem2nFound(i+(j-1)*sum(WallNoElem),6) = 1 + size(wallNodesXYZ,1)...
                   + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
                   + sum(WallNoElem)*(foundationOutNoElem + WallThickNoElem);
           end
       end  
    end
end

for i = 2:foundationHeighNoElem
    elem2nFound(sum(WallNoElem)*(foundationOutNoElem+WallThickNoElem+foundationInNoElem)*(i-1)+1:...
        sum(WallNoElem)*(foundationOutNoElem+WallThickNoElem+foundationInNoElem)*i, 1:4) = ...    
    elem2nFound(sum(WallNoElem)*(foundationOutNoElem+WallThickNoElem+foundationInNoElem)*(i-2)+1:...
        sum(WallNoElem)*(foundationOutNoElem+WallThickNoElem+foundationInNoElem)*(i-1), 5:8);
    
    elem2nFound(sum(WallNoElem)*(foundationOutNoElem+WallThickNoElem+foundationInNoElem)*(i-1)+1:...
        sum(WallNoElem)*(foundationOutNoElem+WallThickNoElem+foundationInNoElem)*i, 5:8) = ...
    elem2nFound(sum(WallNoElem)*(foundationOutNoElem+WallThickNoElem+foundationInNoElem)*(i-1)+1:...
        sum(WallNoElem)*(foundationOutNoElem+WallThickNoElem+foundationInNoElem)*i, 1:4)...
        + sum(WallNoElem)*(foundationOutNoElem+WallThickNoElem+foundationInNoElem+1);    
end

elem2nFoundReorder = elem2nFound;
elem2nFoundReorder(:,1:4) = elem2nFound(:,5:8);
elem2nFoundReorder(:,5:8) = elem2nFound(:,1:4);
wholeElem2n = [elem2n; elem2nFoundReorder];
PlotMesh(wholeNodesXYZ, wholeElem2n,0)
endFoundationElem = length(wholeElem2n);
%% Add internal wall
%---------------------inter wall 1-----------------------------------------
interWall = [1 2];
interOpeningNodes = [8-partWallThick/2, 8+partWallThick/2, ...
    16-partWallThick/2, 16+partWallThick/2, ... 
    24-partWallThick/2, 24+partWallThick/2, ...
    32-partWallThick/2, 32+partWallThick/2, ...
    40-partWallThick/2, 40+partWallThick/2];
interOpeningNodes = interOpeningNodes - WallThick(4);
interOpening = [];
for i  = 1:length(interOpeningNodes)
    interOpening = [interOpening; [1 interOpeningNodes(i) 0 interOpeningNodes(i) 0]];
end
InterNodesXYZ = [];

interWallCorners = [-24+0.2 0-partWallThick/2 0;
                    24-0.2  0-partWallThick/2 0];
[InterCentralNodesXYZ, InterWallNoElem, wallSegNo, segElemNo]= meshOnCurveWithFixNodes(interWallCorners, interWall, ...
    WallExterElemApproSize, interOpening);
InterCentralNodesXYZ(end,:) = [];
InterNodesXYZ = [InterNodesXYZ; InterCentralNodesXYZ];

interWallCorners = [-24+0.2 0+partWallThick/2 0;
                    24-0.2  0+partWallThick/2 0];
[InterCentralNodesXYZ, InterWallNoElem, wallSegNo, segElemNo]= meshOnCurveWithFixNodes(interWallCorners, interWall, ...
    WallExterElemApproSize, interOpening);
InterCentralNodesXYZ(end,:) = [];
InterNodesXYZ = [InterNodesXYZ; InterCentralNodesXYZ];
temp = InterNodesXYZ;
for i = 2:length(heightZ)
    temp(:,3) = heightZ(i);
    InterNodesXYZ = [InterNodesXYZ; temp];
end
for i = 1:foundationHeighNoElem
    temp(:,3) = 0-i*foundationHeight/foundationHeighNoElem;
    InterNodesXYZ = [InterNodesXYZ; temp];
end 

% Elements of internal wall 1 
interElem2n = zeros(InterWallNoElem*1*(WallHeightNoElem+2),8);
for i = 1:InterWallNoElem
    for j = 1:WallHeightNoElem
       if i==1
          interElem2n(1 + (j-1)*InterWallNoElem*1,1) = ...
              find(abs(wholeNodesXYZ(:,1)- (-24+0.2))<1e-10 & wholeNodesXYZ(:,3)== heightZ(j) & abs(wholeNodesXYZ(:,2)-(0-partWallThick/2))<1e-10);
          interElem2n(1 + (j-1)*InterWallNoElem*1,2) = ...
              i + (j-1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
          interElem2n(1 + (j-1)*InterWallNoElem*1,3) = ...
              i + (j-1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1);
          interElem2n(1 + (j-1)*InterWallNoElem*1,4) = ...
              find(abs(wholeNodesXYZ(:,1)- (-24+0.2))<1e-10 & wholeNodesXYZ(:,3)== heightZ(j) & abs(wholeNodesXYZ(:,2)-(0+partWallThick/2))<1e-10);
          interElem2n(1 + (j-1)*InterWallNoElem*1,5) = ...
              find(abs(wholeNodesXYZ(:,1)- (-24+0.2))<1e-10& wholeNodesXYZ(:,3)== heightZ(j+1) & abs(wholeNodesXYZ(:,2)-(0-partWallThick/2))<1e-10);
          interElem2n(1 + (j-1)*InterWallNoElem*1,6) = ...
              i + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
          interElem2n(1 + (j-1)*InterWallNoElem*1,7) = ...
              i + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1);
          interElem2n(1 + (j-1)*InterWallNoElem*1,8) = ...
               find(abs(wholeNodesXYZ(:,1)- (-24+0.2))<1e-10& wholeNodesXYZ(:,3)== heightZ(j+1) & abs(wholeNodesXYZ(:,2)-(0+partWallThick/2))<1e-10);
       elseif i==InterWallNoElem
          interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,1) = ...
              InterWallNoElem-1 + (j-1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
          interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,2) = ...
              find(abs(wholeNodesXYZ(:,1)- (24-0.2))<1e-10 & wholeNodesXYZ(:,3)== heightZ(j) & abs(wholeNodesXYZ(:,2)-(0-partWallThick/2))<1e-10);
          interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,3) = ...
              find(abs(wholeNodesXYZ(:,1)- (24-0.2))<1e-10 & wholeNodesXYZ(:,3)== heightZ(j) & abs(wholeNodesXYZ(:,2)-(0+partWallThick/2))<1e-10);
          interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,4) = ...
              InterWallNoElem-1 + (j-1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1);
          interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,5) = ...
              InterWallNoElem-1 + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
          interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,6) = ...
              find(abs(wholeNodesXYZ(:,1)- (24-0.2))<1e-10 & wholeNodesXYZ(:,3)== heightZ(j+1) & abs(wholeNodesXYZ(:,2)-(0-partWallThick/2))<1e-10);
          interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,7) = ...
              find(abs(wholeNodesXYZ(:,1)- (24-0.2))<1e-10 & wholeNodesXYZ(:,3)== heightZ(j+1) & abs(wholeNodesXYZ(:,2)-(0+partWallThick/2))<1e-10);
          interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,8) = ...
              InterWallNoElem-1 + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1);
       else
          interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,1) = ...
              i-1 + (j-1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
          interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,2) = ...
              i-1 + (j-1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1)+1;
          interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,3) = ...
              i-1 + (j-1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1)+1+(InterWallNoElem-1);
          interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,4) = ...
              i-1 + (j-1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1)+(InterWallNoElem-1);
          interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,5) = ...
              i-1 + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
          interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,6) = ...
              i-1 + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1)+1;
          interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,7) = ...
              i-1 + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1)+1+(InterWallNoElem-1);
          interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,8) = ...
              i-1 + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1)+(InterWallNoElem-1);           
       end
    end
end
% Add foundation elements to Internal wall 1, note this only works for
% foundationHeighNoElem = 2
for i = 1:InterWallNoElem
     j = WallHeightNoElem;
       if i==1
          leftConnect = -24+WallThick(4)+(foundationWidth(4)-WallThick(4))/2;
          rightConnect = -leftConnect;
          interElem2n(1 + (j)*InterWallNoElem*1,1) = ...
              find(abs(wholeNodesXYZ(:,1)- (leftConnect))<1e-10 & wholeNodesXYZ(:,3)== -foundationHeight/2 & abs(wholeNodesXYZ(:,2)-(0-partWallThick/2))<1e-10);
          interElem2n(1 + (j)*InterWallNoElem*1,2) = ...
              i + (j+1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
          interElem2n(1 + (j)*InterWallNoElem*1,3) = ...
              i + (j+1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1);
          interElem2n(1 + (j)*InterWallNoElem*1,4) = ...
              find(abs(wholeNodesXYZ(:,1)- (leftConnect))<1e-10 & wholeNodesXYZ(:,3)== -foundationHeight/2 & abs(wholeNodesXYZ(:,2)-(0+partWallThick/2))<1e-10);
          interElem2n(1 + (j)*InterWallNoElem*1,5) = ...
              find(abs(wholeNodesXYZ(:,1)- (leftConnect))<1e-10& wholeNodesXYZ(:,3)== 0 & abs(wholeNodesXYZ(:,2)-(0-partWallThick/2))<1e-10);
          interElem2n(1 + (j)*InterWallNoElem*1,6) = ...
              i + (0)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
          interElem2n(1 + (j)*InterWallNoElem*1,7) = ...
              i + (0)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1);
          interElem2n(1 + (j)*InterWallNoElem*1,8) = ...
               find(abs(wholeNodesXYZ(:,1)- (leftConnect))<1e-10& wholeNodesXYZ(:,3)== 0 & abs(wholeNodesXYZ(:,2)-(0+partWallThick/2))<1e-10);
       elseif i==InterWallNoElem
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,1) = ...
              InterWallNoElem-1 + (j+1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,2) = ...
              find(abs(wholeNodesXYZ(:,1)- (rightConnect))<1e-10 & wholeNodesXYZ(:,3)== -foundationHeight & abs(wholeNodesXYZ(:,2)-(0-partWallThick/2))<1e-10);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,3) = ...
              find(abs(wholeNodesXYZ(:,1)- (rightConnect))<1e-10 & wholeNodesXYZ(:,3)== -foundationHeight & abs(wholeNodesXYZ(:,2)-(0+partWallThick/2))<1e-10);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,4) = ...
              InterWallNoElem-1 + (j+1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,5) = ...
              InterWallNoElem-1 + (0)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,6) = ...
              find(abs(wholeNodesXYZ(:,1)- (rightConnect))<1e-10 & wholeNodesXYZ(:,3)== 0 & abs(wholeNodesXYZ(:,2)-(0-partWallThick/2))<1e-10);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,7) = ...
              find(abs(wholeNodesXYZ(:,1)- (rightConnect))<1e-10 & wholeNodesXYZ(:,3)== 0 & abs(wholeNodesXYZ(:,2)-(0+partWallThick/2))<1e-10);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,8) = ...
              InterWallNoElem-1 + (0)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1);
       else
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,1) = ...
              i-1 + (j+1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,2) = ...
              i-1 + (j+1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1)+1;
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,3) = ...
              i-1 + (j+1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1)+1+(InterWallNoElem-1);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,4) = ...
              i-1 + (j+1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1)+(InterWallNoElem-1);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,5) = ...
              i-1 + (0)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,6) = ...
              i-1 + (0)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1)+1;
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,7) = ...
              i-1 + (0)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1)+1+(InterWallNoElem-1);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,8) = ...
              i-1 + (0)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1)+(InterWallNoElem-1);           
       end
end

for i = 1:InterWallNoElem
     j = WallHeightNoElem+1;
       if i==1
          leftConnect = -24+WallThick(4)+(foundationWidth(4)-WallThick(4))/2;
          rightConnect = -leftConnect;
          interElem2n(1 + (j)*InterWallNoElem*1,1) = ...
              find(abs(wholeNodesXYZ(:,1)- (leftConnect))<1e-10 & wholeNodesXYZ(:,3)== -foundationHeight & abs(wholeNodesXYZ(:,2)-(0-partWallThick/2))<1e-10);
          interElem2n(1 + (j)*InterWallNoElem*1,2) = ...
              i + (j+1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
          interElem2n(1 + (j)*InterWallNoElem*1,3) = ...
              i + (j+1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1);
          interElem2n(1 + (j)*InterWallNoElem*1,4) = ...
              find(abs(wholeNodesXYZ(:,1)- (leftConnect))<1e-10 & wholeNodesXYZ(:,3)== -foundationHeight & abs(wholeNodesXYZ(:,2)-(0+partWallThick/2))<1e-10);
          interElem2n(1 + (j)*InterWallNoElem*1,5) = ...
              find(abs(wholeNodesXYZ(:,1)- (leftConnect))<1e-10& wholeNodesXYZ(:,3)== -foundationHeight/2 & abs(wholeNodesXYZ(:,2)-(0-partWallThick/2))<1e-10);
          interElem2n(1 + (j)*InterWallNoElem*1,6) = ...
              i + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
          interElem2n(1 + (j)*InterWallNoElem*1,7) = ...
              i + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1);
          interElem2n(1 + (j)*InterWallNoElem*1,8) = ...
               find(abs(wholeNodesXYZ(:,1)- (leftConnect))<1e-10& wholeNodesXYZ(:,3)== -foundationHeight/2 & abs(wholeNodesXYZ(:,2)-(0+partWallThick/2))<1e-10);
       elseif i==InterWallNoElem
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,1) = ...
              InterWallNoElem-1 + (j+1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,2) = ...
              find(abs(wholeNodesXYZ(:,1)- (rightConnect))<1e-10 & wholeNodesXYZ(:,3)== -foundationHeight & abs(wholeNodesXYZ(:,2)-(0-partWallThick/2))<1e-10);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,3) = ...
              find(abs(wholeNodesXYZ(:,1)- (rightConnect))<1e-10 & wholeNodesXYZ(:,3)== -foundationHeight & abs(wholeNodesXYZ(:,2)-(0+partWallThick/2))<1e-10);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,4) = ...
              InterWallNoElem-1 + (j+1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,5) = ...
              InterWallNoElem-1 + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,6) = ...
              find(abs(wholeNodesXYZ(:,1)- (rightConnect))<1e-10 & wholeNodesXYZ(:,3)== -foundationHeight/2 & abs(wholeNodesXYZ(:,2)-(0-partWallThick/2))<1e-10);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,7) = ...
              find(abs(wholeNodesXYZ(:,1)- (rightConnect))<1e-10 & wholeNodesXYZ(:,3)== -foundationHeight/2 & abs(wholeNodesXYZ(:,2)-(0+partWallThick/2))<1e-10);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,8) = ...
              InterWallNoElem-1 + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1);
       else
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,1) = ...
              i-1 + (j+1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,2) = ...
              i-1 + (j+1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1)+1;
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,3) = ...
              i-1 + (j+1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1)+1+(InterWallNoElem-1);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,4) = ...
              i-1 + (j+1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1)+(InterWallNoElem-1);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,5) = ...
              i-1 + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,6) = ...
              i-1 + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1)+1;
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,7) = ...
              i-1 + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1)+1+(InterWallNoElem-1);
          interElem2n(1+(i-1) + (j)*InterWallNoElem*1,8) = ...
              i-1 + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1)+(InterWallNoElem-1);           
       end
end
    

% 
wholeNodesXYZ = [wholeNodesXYZ; InterNodesXYZ];
PlotMesh(wholeNodesXYZ, interElem2n, 0)
wholeElem2n = [wholeElem2n; interElem2n];
PlotMesh(wholeNodesXYZ, wholeElem2n, 0)
% ---------------------------inter wall 2 ---------------------------------
interWallCornersListX = [-16, -8, 0, 8, 16];
interWallCornersListYs = [-8+0.2, 0+partWallThick/2];
interWallCornersListYt = [0-partWallThick/2, 8-0.2];
for interI = 1:length(interWallCornersListX)
   for interJ = 1:length(interWallCornersListYs)
        interWall = [1 2];
        interOpening = [];
        InterNodesXYZ = [];
        interWallCorners = [interWallCornersListX(interI)-partWallThick/2  interWallCornersListYs(interJ), 0;
                            interWallCornersListX(interI)-partWallThick/2  interWallCornersListYt(interJ), 0];
        [InterCentralNodesXYZ, InterWallNoElem, wallSegNo, segElemNo]= meshOnCurveWithFixNodes(interWallCorners, interWall, ...
            WallExterElemApproSize, interOpening);
        InterCentralNodesXYZ(end,:) = [];
        InterNodesXYZ = [InterNodesXYZ; InterCentralNodesXYZ];
        interWallCorners = [interWallCornersListX(interI)+partWallThick/2  interWallCornersListYs(interJ), 0;
                            interWallCornersListX(interI)+partWallThick/2  interWallCornersListYt(interJ), 0];
        [InterCentralNodesXYZ, InterWallNoElem, wallSegNo, segElemNo]= meshOnCurveWithFixNodes(interWallCorners, interWall, ...
            WallExterElemApproSize, interOpening);
        InterCentralNodesXYZ(end,:) = [];
        InterNodesXYZ = [InterNodesXYZ; InterCentralNodesXYZ];
        temp = InterNodesXYZ;
        for i = 2:length(heightZ)
            temp(:,3) = heightZ(i);
            InterNodesXYZ = [InterNodesXYZ; temp];
        end
        interElem2n = zeros(InterWallNoElem*1*(WallHeightNoElem),8);
        for i = 1:InterWallNoElem
            for j = 1:WallHeightNoElem
               if i==1
                  interElem2n(1 + (j-1)*InterWallNoElem*1,1) = ...
                      find(abs(wholeNodesXYZ(:,1)- (interWallCornersListX(interI)-partWallThick/2))<1e-10 & wholeNodesXYZ(:,3)== heightZ(j) & abs(wholeNodesXYZ(:,2)-interWallCornersListYs(interJ))<1e-10);
                  interElem2n(1 + (j-1)*InterWallNoElem*1,2) = ...
                      find(abs(wholeNodesXYZ(:,1)- (interWallCornersListX(interI)+partWallThick/2))<1e-10 & wholeNodesXYZ(:,3)== heightZ(j) & abs(wholeNodesXYZ(:,2)-interWallCornersListYs(interJ))<1e-10);
                  interElem2n(1 + (j-1)*InterWallNoElem*1,3) = ...
                      i + (j-1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1);
                  interElem2n(1 + (j-1)*InterWallNoElem*1,4) = ...
                      i + (j-1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
                  interElem2n(1 + (j-1)*InterWallNoElem*1,5) = ...
                      find(abs(wholeNodesXYZ(:,1)- (interWallCornersListX(interI)-partWallThick/2))<1e-10 & wholeNodesXYZ(:,3)== heightZ(j+1) & abs(wholeNodesXYZ(:,2)-interWallCornersListYs(interJ))<1e-10);
                  interElem2n(1 + (j-1)*InterWallNoElem*1,6) = ...
                      find(abs(wholeNodesXYZ(:,1)- (interWallCornersListX(interI)+partWallThick/2))<1e-10 & wholeNodesXYZ(:,3)== heightZ(j+1) & abs(wholeNodesXYZ(:,2)-interWallCornersListYs(interJ))<1e-10);
                  interElem2n(1 + (j-1)*InterWallNoElem*1,7) = ...
                      i + (j-1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1)*3;
                  interElem2n(1 + (j-1)*InterWallNoElem*1,8) = ...
                      i + (j-1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1)*2;
               elseif i==InterWallNoElem
                  interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,1) = ...
                      (j-1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1);
                  interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,2) = ...
                      (j-1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1)*2;
                  interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,3) = ...
                      find(abs(wholeNodesXYZ(:,1)- (interWallCornersListX(interI)+partWallThick/2))<1e-10 & wholeNodesXYZ(:,3)== heightZ(j) & abs(wholeNodesXYZ(:,2)-interWallCornersListYt(interJ))<1e-10);
                  interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,4) = ...
                      find(abs(wholeNodesXYZ(:,1)- (interWallCornersListX(interI)-partWallThick/2))<1e-10 & wholeNodesXYZ(:,3)== heightZ(j) & abs(wholeNodesXYZ(:,2)-interWallCornersListYt(interJ))<1e-10);
                  interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,5) = ...
                      (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1);
                  interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,6) = ...
                      (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1)*2;
                  interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,7) = ...
                      find(abs(wholeNodesXYZ(:,1)- (interWallCornersListX(interI)+partWallThick/2))<1e-10 & wholeNodesXYZ(:,3)== heightZ(j+1) & abs(wholeNodesXYZ(:,2)-interWallCornersListYt(interJ))<1e-10);
                  interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,8) = ...
                      find(abs(wholeNodesXYZ(:,1)- (interWallCornersListX(interI)-partWallThick/2))<1e-10 & wholeNodesXYZ(:,3)== heightZ(j+1) & abs(wholeNodesXYZ(:,2)-interWallCornersListYt(interJ))<1e-10);
               else
                  interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,1) = ...
                      i-1 + (j-1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
                  interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,2) = ...
                      i-1 + (j-1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1);
                  interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,3) = ...
                      i-1 + (j-1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + 1 + (InterWallNoElem-1);
                  interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,4) = ...
                      i-1 + (j-1)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + 1;
                  interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,5) = ...
                      i-1 + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1);
                  interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,6) = ...
                      i-1 + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + (InterWallNoElem-1);
                  interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,7) = ...
                      i-1 + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + 1 + (InterWallNoElem-1);
                  interElem2n(1+(i-1) + (j-1)*InterWallNoElem*1,8) = ...
                      i-1 + (j)*2*(InterWallNoElem-1) + size(wholeNodesXYZ,1) + 1;           
               end
            end
        end
        wholeNodesXYZ = [wholeNodesXYZ; InterNodesXYZ];
        wholeElem2n = [wholeElem2n; interElem2n];
   end
end
PlotMesh(wholeNodesXYZ, wholeElem2n, 0)

%% InterFace elems
% elem2nInter = elem2nFoundReorder(end-(WallThickNoElem+foundationOutNoElem+foundationInNoElem)*...
%     sum(WallNoElem)+1:end,:);
% PlotMesh(wholeNodesXYZ, elem2nInter, 0)
wholeNodesXYZ(:,3) = wholeNodesXYZ(:,3)+foundationHeight;
bottom_nodes = find(wholeNodesXYZ(:,3) == min(wholeNodesXYZ(:,3)));

% figure
% plot(wholeNodesXYZ(bottom_nodes,1), wholeNodesXYZ(bottom_nodes,2),'.')
bottom_nodesCenterWall = bottom_nodes(1061:end);
bottom_nodes = bottom_nodes(1:1060);
bottomNodesXYZ = wholeNodesXYZ(bottom_nodes,:);
outerDummyNodesXYZ = offsetWallNodes(wallNodesXYZ(1:sum(WallNoElem),:),...
    wallNodesXYZ(1+sum(WallNoElem):sum(WallNoElem)*2,:),WallThick(1)/WallThickNoElem,...
    (foundationWidth(1)-WallThick(1))/2 +(foundationWidth(1)-WallThick(1))/2/foundationOutNoElem);
innerDummyNodesXYZ = offsetWallNodes(wallNodesXYZ(1+sum(WallNoElem)*WallThickNoElem:sum(WallNoElem)+sum(WallNoElem)*WallThickNoElem,:),...
        wallNodesXYZ(1+sum(WallNoElem)*(WallThickNoElem-1):sum(WallNoElem)*WallThickNoElem,:),WallThick(1)/WallThickNoElem,...
        (foundationWidth(1)-WallThick(1))/2+(foundationWidth(1)-WallThick(1))/2/foundationInNoElem);
interNodesXYZ = [outerDummyNodesXYZ;bottomNodesXYZ;innerDummyNodesXYZ];

centerWalllower = wholeNodesXYZ(bottom_nodesCenterWall(1:52),:);
centerWallupper = wholeNodesXYZ(bottom_nodesCenterWall(53:end),:);
centerWalllower = [[-24, -0.02, 0];centerWalllower;[24, -0.02, 0]]; 
centerWallupper = [[-24, 0.02, 0];centerWallupper;[24, 0.02, 0]]; 
centerWallXYZ = [centerWalllower; centerWallupper];
centerWalllower(:,2) = centerWalllower(:,2) - 0.04;
centerWallupper(:,2) = centerWallupper(:,2) + 0.04;
interNodesXYZ = [interNodesXYZ;centerWalllower;centerWallXYZ;centerWallupper];

elem2nInter = zeros((WallThickNoElem+foundationOutNoElem+foundationInNoElem+2)*...
    sum(WallNoElem),4);
for j= 1:foundationOutNoElem+WallThickNoElem+foundationInNoElem+2
    for i = 1:sum(WallNoElem)
       elem2nInter(i+(j-1)*sum(WallNoElem),1) = i+(j-1)*sum(WallNoElem);
       elem2nInter(i+(j-1)*sum(WallNoElem),2) = i+1 + (j-1)*sum(WallNoElem);
       elem2nInter(i+(j-1)*sum(WallNoElem),3) = i+1 + sum(WallNoElem) + (j-1)*sum(WallNoElem);
       elem2nInter(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem);
       if i==sum(WallNoElem)
           elem2nInter(i+(j-1)*sum(WallNoElem),1) = i + (j-1)*sum(WallNoElem);
           elem2nInter(i+(j-1)*sum(WallNoElem),2) = 1 + (j-1)*sum(WallNoElem);
           elem2nInter(i+(j-1)*sum(WallNoElem),3) = 1 + sum(WallNoElem) + (j-1)*sum(WallNoElem);
           elem2nInter(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem);
       end  
    end
end
% PlotMesh(interNodesXYZ(:,1:2), elem2nInter, 0)
elem2nInterCenterWall = zeros(53*3,4);
for i=1:53
    for j=0:2
        elem2nInterCenterWall(i+j*53,1) = length(interNodesXYZ)-54*4 + i +j*54;
        elem2nInterCenterWall(i+j*53,2) = length(interNodesXYZ)-54*4 + i+1 +j*54;
        elem2nInterCenterWall(i+j*53,3) = length(interNodesXYZ)-54*4 + i+1+54 +j*54;
        elem2nInterCenterWall(i+j*53,4) = length(interNodesXYZ)-54*4 + i+54+j*54;
    end
end
% PlotMesh(interNodesXYZ(:,1:2), elem2nInterCenterWall, 0) 
elem2nInter = [elem2nInter;elem2nInterCenterWall];
PlotMesh(interNodesXYZ(:,1:2), elem2nInter, 0)
%% Coordinate transfer:
% move perpendicular to y axis by 24m
wholeNodesXYZ(:,1) = wholeNodesXYZ(:,1) - 24;
interNodesXYZ(:,1) = interNodesXYZ(:,1) - 24;

% move along y axis by -8m
wholeNodesXYZ(:,2) = wholeNodesXYZ(:,2) + 8;
interNodesXYZ(:,2) = interNodesXYZ(:,2) + 8;

% rotate clockwise for 26 degree
theta = 26/180*pi;
R = [cos(theta) sin(theta); -sin(theta), cos(theta)];
for i = 1:length(wholeNodesXYZ)
    wholeNodesXYZ(i, 1:2) = (R*wholeNodesXYZ(i, 1:2)')';
end
for i = 1:length(interNodesXYZ)
    interNodesXYZ(i, 1:2) = (R*interNodesXYZ(i, 1:2)')';
end
PlotMesh(wholeNodesXYZ, wholeElem2n,0)
PlotMesh(interNodesXYZ(:,1:2), elem2nInter, 0)
PlotMesh(wholeNodesXYZ, wholeElem2n(startFoundationElem:endFoundationElem,:),0)
%% Save
% oldFolder = cd(oldFolder);
timeString = datestr(datetime('now'),30);
timberElemIndex = [1 1];
title = strcat('camos_mat_',timeString,'.mat');
save(title,'wholeNodesXYZ','wholeElem2n','interNodesXYZ','elem2nInter', 'timberElemIndex')
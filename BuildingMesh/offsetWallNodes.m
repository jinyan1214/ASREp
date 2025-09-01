% function newNodesXYZ = offsetWallNodes(nodes1, nodes2, basedist, offdist,...
%     WallNoElem)
% %     distDiff = sqrt(sum((nodes1(:,1:2)-nodes2(:,1:2)).^2,2));
% %     newx = nodes1(:,1)+(nodes1(:,1)-nodes2(:,1))/basedist*offdist;
% %     newy = nodes1(:,2)+(nodes1(:,2)-nodes2(:,2))/basedist*offdist;
% %     newNodesXYZ = [newx, newy, nodes1(:,3)];
%     
%     newx = zeros(size(nodes1(:,1)));
%     newy = zeros(size(nodes1(:,1)));
%     for i = 1:length(WallNoElem)
%         startInd = sum(WallNoElem(1:i-1))+1;
%         endInd = sum(WallNoElem(1:i));
%         newx(startInd:endInd) = nodes1(startInd:endInd,1)+...
%             (nodes1(startInd:endInd,1)-nodes2(startInd:endInd,1))/...
%             basedist(i)*offdist(i);
%         newy(startInd:endInd) = nodes1(startInd:endInd,2)+...
%             (nodes1(startInd:endInd,2)-nodes2(startInd:endInd,2))/...
%             basedist(i)*offdist(i);
%     end
%     newNodesXYZ = [newx, newy, nodes1(:,3)];
% end
function newNodesXYZ = offsetWallNodes(nodes1, nodes2, basedist, offdist)
%     distDiff = sqrt(sum((nodes1(:,1:2)-nodes2(:,1:2)).^2,2));
    newx = nodes1(:,1)+(nodes1(:,1)-nodes2(:,1))/basedist*offdist;
    newy = nodes1(:,2)+(nodes1(:,2)-nodes2(:,2))/basedist*offdist;
    newNodesXYZ = [newx, newy, nodes1(:,3)];
end

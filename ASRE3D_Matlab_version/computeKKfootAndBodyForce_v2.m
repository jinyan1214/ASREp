function [KKfoot, Q] = computeKKfootAndBodyForce_v2(nX, nY, nZ, elem2n, ...
    elemE, elemv, elemRho, elemG)
% nX: x coordinate of nodes
% nY: y coordinate of nodes
% nZ: z coordinate of nodes
% elem2n: n_elem by 8 matrix(i,j), # of jth node in the ith element, the
% nodes are ordered counter-clockwise, bottom to up
% elemE: young's modulus of the elements
% elemv: poission's ratio of the elements
% The nth node, the global DOF: 3n-2, 3n-1, 3n

dimKK = 3 * size(nX, 1);
if dimKK>50000
    KKfoot = sparse(dimKK, dimKK);
else
    KKfoot = zeros(dimKK, dimKK);
end
Q = zeros(dimKK, 1);
for i = 1:size(elem2n, 1)
   xyzLocal = zeros(8, 3);
   nGlobal = [];
   nGlobalVertical = [];
    for j = 1:8
        xyzLocal(j,:) = [nX(elem2n(i,j)), nY(elem2n(i,j)), nZ(elem2n(i,j))];
        nGlobal = [nGlobal,3*elem2n(i,j)-2,3*elem2n(i,j)-1,3*elem2n(i,j)];
        nGlobalVertical = [nGlobalVertical,3*elem2n(i,j)];
    end
    [kkLocal, vol] = computeKLocalLE8nodeBrickAndVolume_v2(xyzLocal, elemE(i), elemv(i), elemG(i));
    Q(nGlobalVertical) = Q(nGlobalVertical) - elemRho(i) * vol/8;
    KKfoot(nGlobal, nGlobal) = KKfoot(nGlobal, nGlobal) + kkLocal;  
end

return
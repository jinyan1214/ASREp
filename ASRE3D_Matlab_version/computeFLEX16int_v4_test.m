function FLEX = computeFLEX16int_v4_test(XYZ, elem2n, connectXY, Gs, nus)
    realNodes = unique(reshape(elem2n, length(elem2n)*4,1));
    dim = length(realNodes)*3;
    FLEX = zeros(dim, dim);
    FLEX2 = zeros(dim,dim);
    realXYZ = XYZ(realNodes,:);
    realConnect = connectXY(realNodes,:);
    nat1 = sqrt(3/7-2/7*sqrt(6/5));
    nat2 = sqrt(3/7+2/7*sqrt(6/5));
    nat  = [ -nat2 -nat2; -nat1 -nat2; +nat1 -nat2; +nat2 -nat2;...
        -nat2 -nat1; -nat1 -nat1; +nat1 -nat1; +nat2 -nat1;...
        -nat2 +nat1; -nat1 +nat1; +nat1 +nat1; +nat2 +nat1;...
        -nat2 +nat2; -nat1 +nat2; +nat1 +nat2; +nat2 +nat2;...
            ];
    wi1 = (18+sqrt(30))/36;
    wi2 = (18-sqrt(30))/36;
    wi = [wi2*wi2 wi2*wi1 wi2*wi1 wi2*wi2...
        wi1*wi2 wi1*wi1 wi1*wi1 wi1*wi2...
        wi1*wi2 wi1*wi1 wi1*wi1 wi1*wi2...
        wi2*wi2 wi2*wi1 wi2*wi1 wi2*wi2];
    nodeCount = [0, 0];
    inOrOutVec = zeros(length(realNodes),1);
    for j=1:length(realNodes)
        conx = realConnect(j,1:2:end);
        cony = realConnect(j,2:2:end);
        [conxS, conyS] = sortCon(conx, cony);
        if insidePoly(realXYZ(j,1), realXYZ(j,2), conxS, conyS)
            inOrOutVec(j) = 1;
            nodeCount(1) = nodeCount(1) + 1;
            isSouth = cony<realXYZ(j,2);
            isWest = conx<realXYZ(j,1);
%             if sum(isSouth)~=2 || sum(isWest)~=2
%                 patchL = realXYZ(j,1)-0.05/2;
%                 patchR = realXYZ(j,1)+0.05/2;
%                 patchB = realXYZ(j,2)-0.05/2;
%                 patchT = realXYZ(j,2)+0.05/2;
%                 [ uxij_1,uyij_1,uzij_1 ] = int_factor_mindlin_j_1_vect_Vasiri_jz_boundary(realXYZ(:,1),realXYZ(:,2),realXYZ(:,3),...
%                     realXYZ(j,1),realXYZ(j,2),realXYZ(j,3),...
%                     Gs*ones(length(realNodes),1),nus*ones(length(realNodes),1),...
%                     patchL, patchR, patchB, patchT);
%                 [ uxij_2,uyij_2,uzij_2 ] = int_factor_mindlin_j_2_vect_Vasiri_jz_boundary(realXYZ(:,1),realXYZ(:,2),realXYZ(:,3),...
%                     realXYZ(j,1),realXYZ(j,2),realXYZ(j,3),...
%                     Gs*ones(length(realNodes),1),nus*ones(length(realNodes),1),...
%                     patchL, patchR, patchB, patchT);
%                 [ uxij_3,uyij_3,uzij_3 ] = int_factor_mindlin_j_3_vect_Vasiri_boundary(realXYZ(:,1),realXYZ(:,2),realXYZ(:,3),...
%                     realXYZ(j,1),realXYZ(j,2),realXYZ(j,3),...
%                     Gs*ones(length(realNodes),1),nus*ones(length(realNodes),1),...
%                     patchL, patchR, patchB, patchT);
%                 if sum(isnan(uxij_1))>0 || sum(isnan(uyij_2))>0 || sum(isnan(uzij_3))>0
%                     numNan = sum(isnan(uxij_1));
%                    uxij_1 = fillmissing(uxij_1, 'movmean',numNan+1);
%                    uxij_2 = fillmissing(uxij_2, 'movmean',numNan+1);
%                    uxij_3 = fillmissing(uxij_3, 'movmean',numNan+1);
%                    uyij_1 = fillmissing(uyij_1, 'movmean',numNan+1);
%                    uyij_2 = fillmissing(uyij_2, 'movmean',numNan+1);
%                    uyij_3 = fillmissing(uyij_3, 'movmean',numNan+1);
%                    uzij_1 = fillmissing(uzij_1, 'movmean',numNan+1);
%                    uzij_2 = fillmissing(uzij_2, 'movmean',numNan+1);
%                    uzij_3 = fillmissing(uzij_3, 'movmean',numNan+1);
%                 end
%                 FLEX(1:3:end, j*3-2) = uxij_1;FLEX(1:3:end, j*3-1) = uyij_1;FLEX(1:3:end, j*3-0) = uzij_1;
%                 FLEX(2:3:end, j*3-2) = uxij_2;FLEX(2:3:end, j*3-1) = uyij_2;FLEX(2:3:end, j*3-0) = uzij_2;
%                 FLEX(3:3:end, j*3-2) = uxij_3;FLEX(3:3:end, j*3-1) = uyij_3;FLEX(3:3:end, j*3-0) = uzij_3;
% 
%                 FLEX2(1:3:end, j*3-2) = uxij_1;FLEX2(1:3:end, j*3-1) = uyij_1;FLEX2(1:3:end, j*3-0) = uzij_1;
%                 FLEX2(2:3:end, j*3-2) = uxij_2;FLEX2(2:3:end, j*3-1) = uyij_2;FLEX2(2:3:end, j*3-0) = uzij_2;
%                 FLEX2(3:3:end, j*3-2) = uxij_3;FLEX2(3:3:end, j*3-1) = uyij_3;FLEX2(3:3:end, j*3-0) = uzij_3;
%             else        
%                 sortedx(1) = conx(isSouth&isWest);
%                 sortedy(1) = cony(isSouth&isWest);
%                 sortedx(2) = conx(isSouth&(~isWest));
%                 sortedy(2) = cony(isSouth&(~isWest));
%                 sortedx(3) = conx((~isSouth)&(~isWest));
%                 sortedy(3) = cony((~isSouth)&(~isWest));
%                 sortedx(4) = conx((~isSouth)&(isWest));
%                 sortedy(4) = cony((~isSouth)&(isWest));
                
                sortedx(1) = conxS(1);
                sortedy(1) = conyS(1);
                sortedx(2) = conxS(2);
                sortedy(2) = conyS(2);
                sortedx(3) = conxS(3);
                sortedy(3) = conyS(3);
                sortedx(4) = conxS(4);
                sortedy(4) = conyS(4);

                vecA = [sortedx(3) - sortedx(1), sortedy(3) - sortedy(1), 0];
                vecB = [sortedx(4) - sortedx(2), sortedy(4) - sortedy(2), 0];
    %         patch_A = norm(cross(vecA, vecB))/2/(4*patch_A/4); %area A/4
                patch_area = norm(cross(vecA, vecB))/2;
                patch_A = 1/4; %area A/4
                duijx_1 = 0; duijx_2 = 0; duijx_3 = 0;
                duijy_1 = 0; duijy_2 = 0; duijy_3 = 0;
                duijz_1 = 0; duijz_2 = 0; duijz_3 = 0;
                    for int_i = 1:16
                       xi = nat(int_i,1);
                       eta = nat(int_i,2);
                       N(1) = 1/4 * (1-xi)*(1-eta);
                       N(2) = 1/4 * (1+xi)*(1-eta);
                       N(3) = 1/4 * (1+xi)*(1+eta);
                       N(4) = 1/4 * (1-xi)*(1+eta);
                       x_int = 0; y_int = 0;
                       for k=1:4
                          x_int =  x_int + N(k) * sortedx(k);
                          y_int =  y_int + N(k) * sortedy(k);
                       end
        %                 plot(x_int, y_int, 'rx')
                        [ uxij_1,uyij_1,uzij_1 ] = point_factor_mindlin_j_1_vect_Vasiri_jz(realXYZ(:,1),realXYZ(:,2),realXYZ(:,3),...
                            x_int,y_int,zeros(length(realNodes),1),Gs(j),nus(j));

                        [ uxij_2,uyij_2,uzij_2 ] = point_factor_mindlin_j_2_vect_Vasiri_jz(realXYZ(:,1),realXYZ(:,2),realXYZ(:,3),...
                            x_int,y_int,zeros(length(realNodes),1),Gs(j),nus(j));

                        [ uxij_3,uyij_3,uzij_3 ] = point_factor_mindlin_j_3_vect_Vasiri(realXYZ(:,1),realXYZ(:,2),realXYZ(:,3),...
                            x_int,y_int,zeros(length(realNodes),1),Gs(j),nus(j));
%                         [ uxij_1,uyij_1,uzij_1, uxij_2,uyij_2,uzij_2, uxij_3,uyij_3,uzij_3] = ...
%                             uxijFillNan([uxij_1,uyij_1,uzij_1, uxij_2,uyij_2,uzij_2, uxij_3,uyij_3,uzij_3]);
                        duijx_1 = duijx_1 + uxij_1*wi(int_i);duijy_1 = duijy_1 + uyij_1*wi(int_i);duijz_1 = duijz_1 + uzij_1*wi(int_i);
                        duijx_2 = duijx_2 + uxij_2*wi(int_i);duijy_2 = duijy_2 + uyij_2*wi(int_i);duijz_2 = duijz_2 + uzij_2*wi(int_i);
                        duijx_3 = duijx_3 + uxij_3*wi(int_i);duijy_3 = duijy_3 + uyij_3*wi(int_i);duijz_3 = duijz_3 + uzij_3*wi(int_i);
                    end
                FLEX(1:3:end, j*3-2) = duijx_1*patch_A;FLEX(1:3:end, j*3-1) = duijy_1*patch_A;FLEX(1:3:end, j*3-0) = duijz_1*patch_A;
                FLEX(2:3:end, j*3-2) = duijx_2*patch_A;FLEX(2:3:end, j*3-1) = duijy_2*patch_A;FLEX(2:3:end, j*3-0) = duijz_2*patch_A;
                FLEX(3:3:end, j*3-2) = duijx_3*patch_A;FLEX(3:3:end, j*3-1) = duijy_3*patch_A;FLEX(3:3:end, j*3-0) = duijz_3*patch_A;

                patch_l = (norm([sortedx(1),sortedy(1)]-[sortedx(2),sortedy(2)]) + ...
                    norm([sortedx(4),sortedy(4)]-[sortedx(3),sortedy(3)]))/2;
                patch_h = (norm([sortedx(2),sortedy(2)]-[sortedx(3),sortedy(3)]) + ...
                    norm([sortedx(4),sortedy(4)]-[sortedx(1),sortedy(1)]))/2;
                patchL = realXYZ(j,1)-patch_l/2;
                patchR = realXYZ(j,1)+patch_l/2;
                patchB = realXYZ(j,2)-patch_h/2;
                patchT = realXYZ(j,2)+patch_h/2;
                [ uxij_1,uyij_1,uzij_1 ] = int_factor_mindlin_j_1_vect_Vasiri_jz_boundary(realXYZ(j,1),realXYZ(j,2),realXYZ(j,3),...
                    realXYZ(j,1),realXYZ(j,2),realXYZ(j,3),...
                    Gs(j),nus(j),...
                    patchL, patchR, patchB, patchT);
                [ uxij_2,uyij_2,uzij_2 ] = int_factor_mindlin_j_2_vect_Vasiri_jz_boundary(realXYZ(j,1),realXYZ(j,2),realXYZ(j,3),...
                    realXYZ(j,1),realXYZ(j,2),realXYZ(j,3),...
                    Gs(j),nus(j),...
                    patchL, patchR, patchB, patchT);
                [ uxij_3,uyij_3,uzij_3 ] = int_factor_mindlin_j_3_vect_Vasiri_boundary(realXYZ(j,1),realXYZ(j,2),realXYZ(j,3),...
                    realXYZ(j,1),realXYZ(j,2),realXYZ(j,3),...
                    Gs(j),nus(j),...
                    patchL, patchR, patchB, patchT);
                FLEX(j*3-2, j*3-2) = uxij_1;FLEX(j*3-2, j*3-1) = uyij_1;FLEX(j*3-2, j*3-0) = uzij_1;
                FLEX(j*3-1, j*3-2) = uxij_2;FLEX(j*3-1, j*3-1) = uyij_2;FLEX(j*3-1, j*3-0) = uzij_2;
                FLEX(j*3-0, j*3-2) = uxij_3;FLEX(j*3-0, j*3-1) = uyij_3;FLEX(j*3-0, j*3-0) = uzij_3;
%             end
%             patchL = (sortedx(1)+sortedx(4))/2;
%             patchR = (sortedx(2)+sortedx(3))/2;
%             patchB = (sortedy(1)+sortedy(2))/2;
%             patchT = (sortedy(3)+sortedy(4))/2;
%             [ uxij_1,uyij_1,uzij_1 ] = int_factor_mindlin_j_1_vect_Vasiri_jz_boundary(realXYZ(:,1),realXYZ(:,2),realXYZ(:,3),...
%                 realXYZ(j,1),realXYZ(j,2),realXYZ(j,3),...
%                 Gs*ones(length(realNodes),1),nus*ones(length(realNodes),1),...
%                 patchL, patchR, patchB, patchT);
%             [ uxij_2,uyij_2,uzij_2 ] = int_factor_mindlin_j_2_vect_Vasiri_jz_boundary(realXYZ(:,1),realXYZ(:,2),realXYZ(:,3),...
%                 realXYZ(j,1),realXYZ(j,2),realXYZ(j,3),...
%                 Gs*ones(length(realNodes),1),nus*ones(length(realNodes),1),...
%                 patchL, patchR, patchB, patchT);
%             [ uxij_3,uyij_3,uzij_3 ] = int_factor_mindlin_j_3_vect_Vasiri_boundary(realXYZ(:,1),realXYZ(:,2),realXYZ(:,3),...
%                 realXYZ(j,1),realXYZ(j,2),realXYZ(j,3),...
%                 Gs*ones(length(realNodes),1),nus*ones(length(realNodes),1),...
%                 patchL, patchR, patchB, patchT);
%             if sum(isnan(uxij_1))>0 || sum(isnan(uyij_2))>0 || sum(isnan(uzij_3))>0
%                 numNan = sum(isnan(uxij_1));
%                uxij_1 = fillmissing(uxij_1, 'movmean',numNan+1);
%                uxij_2 = fillmissing(uxij_2, 'movmean',numNan+1);
%                uxij_3 = fillmissing(uxij_3, 'movmean',numNan+1);
%                uyij_1 = fillmissing(uyij_1, 'movmean',numNan+1);
%                uyij_2 = fillmissing(uyij_2, 'movmean',numNan+1);
%                uyij_3 = fillmissing(uyij_3, 'movmean',numNan+1);
%                uzij_1 = fillmissing(uzij_1, 'movmean',numNan+1);
%                uzij_2 = fillmissing(uzij_2, 'movmean',numNan+1);
%                uzij_3 = fillmissing(uzij_3, 'movmean',numNan+1);
%             end
%             FLEX2(1:3:end, j*3-2) = uxij_1;FLEX2(1:3:end, j*3-1) = uyij_1;FLEX2(1:3:end, j*3-0) = uzij_1;
%             FLEX2(2:3:end, j*3-2) = uxij_2;FLEX2(2:3:end, j*3-1) = uyij_2;FLEX2(2:3:end, j*3-0) = uzij_2;
%             FLEX2(3:3:end, j*3-2) = uxij_3;FLEX2(3:3:end, j*3-1) = uyij_3;FLEX2(3:3:end, j*3-0) = uzij_3;
        else
            nodeCount(2) = nodeCount(2) + 1;
            patchSize = 0.05;
            patchL = realXYZ(j,1)-patchSize/2;
            patchR = realXYZ(j,1)+patchSize/2;
            patchB = realXYZ(j,2)-patchSize/2;
            patchT = realXYZ(j,2)+patchSize/2;
            [ uxij_1,uyij_1,uzij_1 ] = int_factor_mindlin_j_1_vect_Vasiri_jz_boundary(realXYZ(:,1),realXYZ(:,2),realXYZ(:,3),...
                realXYZ(j,1),realXYZ(j,2),realXYZ(j,3),...
                Gs(j),nus(j),...
                patchL, patchR, patchB, patchT);
            [ uxij_2,uyij_2,uzij_2 ] = int_factor_mindlin_j_2_vect_Vasiri_jz_boundary(realXYZ(:,1),realXYZ(:,2),realXYZ(:,3),...
                realXYZ(j,1),realXYZ(j,2),realXYZ(j,3),...
                Gs(j),nus(j),...
                patchL, patchR, patchB, patchT);
            [ uxij_3,uyij_3,uzij_3 ] = int_factor_mindlin_j_3_vect_Vasiri_boundary(realXYZ(:,1),realXYZ(:,2),realXYZ(:,3),...
                realXYZ(j,1),realXYZ(j,2),realXYZ(j,3),...
                Gs(j),nus(j),...
                patchL, patchR, patchB, patchT);
            if sum(isnan(uxij_1))>0 || sum(isnan(uyij_2))>0 || sum(isnan(uzij_3))>0
                numNan = sum(isnan(uxij_1));
               uxij_1 = fillmissing(uxij_1, 'movmean',numNan+1);
               uxij_2 = fillmissing(uxij_2, 'movmean',numNan+1);
               uxij_3 = fillmissing(uxij_3, 'movmean',numNan+1);
               uyij_1 = fillmissing(uyij_1, 'movmean',numNan+1);
               uyij_2 = fillmissing(uyij_2, 'movmean',numNan+1);
               uyij_3 = fillmissing(uyij_3, 'movmean',numNan+1);
               uzij_1 = fillmissing(uzij_1, 'movmean',numNan+1);
               uzij_2 = fillmissing(uzij_2, 'movmean',numNan+1);
               uzij_3 = fillmissing(uzij_3, 'movmean',numNan+1);
            end
            FLEX(1:3:end, j*3-2) = uxij_1;FLEX(1:3:end, j*3-1) = uyij_1;FLEX(1:3:end, j*3-0) = uzij_1;
            FLEX(2:3:end, j*3-2) = uxij_2;FLEX(2:3:end, j*3-1) = uyij_2;FLEX(2:3:end, j*3-0) = uzij_2;
            FLEX(3:3:end, j*3-2) = uxij_3;FLEX(3:3:end, j*3-1) = uyij_3;FLEX(3:3:end, j*3-0) = uzij_3;

%             FLEX2(1:3:end, j*3-2) = uxij_1;FLEX2(1:3:end, j*3-1) = uyij_1;FLEX2(1:3:end, j*3-0) = uzij_1;
%             FLEX2(2:3:end, j*3-2) = uxij_2;FLEX2(2:3:end, j*3-1) = uyij_2;FLEX2(2:3:end, j*3-0) = uzij_2;
%             FLEX2(3:3:end, j*3-2) = uxij_3;FLEX2(3:3:end, j*3-1) = uyij_3;FLEX2(3:3:end, j*3-0) = uzij_3;
        end
        
    end
%     temp = FLEX;
%     FLEX = FLEX2;
%     for i = 1:size(FLEX,1)/3
%        FLEX(1+(i-1)*3:i*3,1+(i-1)*3:i*3) =  temp(1+(i-1)*3:i*3,1+(i-1)*3:i*3);
%     end

    
end


function [ uxij_1,uyij_1,uzij_1, uxij_2,uyij_2,uzij_2, uxij_3,uyij_3,uzij_3] = ...
                        uxijFillNan(matrix)
   for i = 1:size(matrix,2)
       if sum(isnan(matrix(:,i)))>0
           numNan = sum(isnan(matrix(:,i)));
           matrix(:,i) = fillmissing(matrix(:,i), 'movmean',numNan+1);
       end
   end
   uxij_1 = matrix(:,1);
   uyij_1 = matrix(:,2);
   uzij_1 = matrix(:,3);
   uxij_2 = matrix(:,4);
   uyij_2 = matrix(:,5);
   uzij_2 = matrix(:,6);
   uxij_3 = matrix(:,7);
   uyij_3 = matrix(:,8);
   uzij_3 = matrix(:,9);
end

% https://math.stackexchange.com/questions/978642/how-to-sort-vertices-of-a-polygon-in-counter-clockwise-order
function [conxS, conyS] = sortCon(conx, cony)
    centerX = mean(conx);
    centerY = mean(cony);
    theta_vec = zeros(length(conx),1);
    for i = 1:length(conx)
        if cony(i)>centerY
           v1 = [1, 0];
           v2 = [conx(i)-centerX, cony(i)-centerY];
           theta = acos((v1*v2')/(norm(v1)*norm(v2)));
        else
            v1 = [1, 0];
            v2 = [conx(i)-centerX, cony(i)-centerY];
            theta = 2*pi - acos((v1*v2')/(norm(v1)*norm(v2)));
        end
        theta_vec(i) = theta;
    end
    XYtheta = [conx', cony', theta_vec];
    XYtheta = sortrows(XYtheta, 3);
    conxS = (XYtheta(:,1))';
    conyS = (XYtheta(:,2))';
end

% https://stackoverflow.com/questions/1119627/how-to-test-if-a-point-is-inside-of-a-convex-polygon-in-2d-integer-coordinates
function isIn = insidePoly(nx, ny, cx, cy)
    previous_side = nan;
    n_vertices = length(cx);
    for n = 1:n_vertices
       a = [cx(n), cy(n)];
       n_next = n+1;
       if n_next > n_vertices
           n_next = 1;
       end
       b = [cx(n_next), cy(n_next)];
       affine_segment = [b(1)-a(1), b(2)-a(2)];
       affine_point = [nx-a(1), ny-a(2)];
       cosSign = sign(affine_segment(1)*affine_point(2) - affine_segment(2)*affine_point(1));
       if cosSign == 0
           isIn = 0;
           return
       elseif isnan(previous_side)
           previous_side = cosSign;
       elseif previous_side ~= cosSign
           isIn = 0;
           return
       end
    end
    isIn = 1;
end
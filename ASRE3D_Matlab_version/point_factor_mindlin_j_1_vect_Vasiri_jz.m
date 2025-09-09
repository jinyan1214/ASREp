function [ uxij_1,uyij_1,uzij_1 ] = point_factor_mindlin_j_1_vect_Vasiri_jz( Xi,Yi,Zi,Xj,Yj,Zj,Gs,nus) 
% uxij_3,uyij_3,uzij_3 displacement at point i due to unit force in j
% in direction 3 (vertical one)
%Xi,Yi,Zi,Xj,Yj,Zj global coordinates of i and j point

%To account that Vaziri and Sun y axis is opposite to Midlin and my RF
Yi=-Yi;
Yj=-Yj;

%Soil properties
ni=nus;
Es=Gs.*2.*(1+ni);
beta=8.*pi.*Es.*(1-ni)./(1+ni);

%Coordinates of the pressure distribution location

X = Xi - Xj;
Y = Yi - Yj;
Z1 = Zi - Zj;
Z2 = Zi + Zj;
R1 = sqrt(X.^2 + Y.^2 + Z1.^2);
R2 = sqrt(X.^2 + Y.^2 + Z2.^2);



delta_x = 1./beta.*((3-4*ni).*(1./R1 + X.^2./(R2.^3)) + 1./R2 + X.^2./(R1.^3) + ...
    2*Zj.*Zi./(R2.^3).*(1 - 3*X.^2./(R2.^2)) + (4*(1-ni).*(1-2*ni)./(R2 + Z2)).*(1 - X.^2./R2./(R2 + Z2)));

delta_y = X.*Y./beta.*(1./(R1.^3) + (3-4*ni)./(R2.^3) - 6.*Zj.*Zi./(R2.^5) - 4*(1-ni).*(1-2*ni)./R2./(R2 + Z2).^2);

delta_z = X./beta.*((3-4*ni).*Z1./(R2.^3) + Z1./(R1.^3) - 6*Zj.*Zi.*Z2./(R2.^5) + 4*(1-ni).*(1-2*ni)./R2./(R2 + Z2));

uxij_1=delta_x;
uyij_1=-delta_y;
uzij_1=delta_z;





end


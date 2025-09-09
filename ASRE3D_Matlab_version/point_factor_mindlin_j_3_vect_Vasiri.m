function [ uxij_3,uyij_3,uzij_3 ] = point_factor_mindlin_j_3_vect_Vasiri( Xi,Yi,Zi,Xj,Yj,Zj,Gs,nus) 
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

%Coordinates of the poit at which the displacement due to the pressure is
%computed
x=Xi;
y=Yi;
z=Zi;

u = Xj;
v = Yj;
w = Zj;

X = x - u;
Y = y - v;
Z1 = z - w;
Z2 = z + w;
R1 = sqrt(X.^2 + Y.^2 + Z1.^2);
R2 = sqrt(X.^2 + Y.^2 + Z2.^2);

delta_z = 1./beta.*((3-4*ni).*(1./R1 + Z2.^2./(R2.^3)) - 2 * w.*z./(R2.^3) + ...
    Z1.^2./(R1.^3) + 6*w.*z.*Z2.^2./(R2.^5) + (8*(1-ni).^2 - (3-4*ni))./R2);
delta_x = X./beta.*((3-4*ni).*Z1./(R2.^3) + Z1./(R1.^3) + 6*w.*z.*Z2./(R2.^5) - ...
    4*(1-ni).*(1-2*ni)./R2./(R2 + Z2));

uzij_3=delta_z;
uxij_3=delta_x;

x=Yi;
y=Xi;
z=Zi;

u = Yj;
v = Xj;
w = Zj;

X = x - u;
Y = y - v;
Z1 = z - w;
Z2 = z + w;
R1 = sqrt(X.^2 + Y.^2 + Z1.^2);
R2 = sqrt(X.^2 + Y.^2 + Z2.^2);

delta_y = X./beta.*((3-4*ni).*Z1./(R2.^3) + Z1./(R1.^3) + 6*w.*z.*Z2./(R2.^5) - ...
    4*(1-ni).*(1-2*ni)./R2./(R2 + Z2));

uyij_3=-delta_y;
%-------------------------------------------------------------------



% uyij_3=-1*0.*y;







end


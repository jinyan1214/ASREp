function [ uxij_3,uyij_3,uzij_3 ] = int_factor_mindlin_j_3_vect_Vasiri_boundary( Xi,Yi,Zi,Xj,Yj,Zj,Gs,nus, patchL, patchR, patchB, patchT) 
% uxij_3,uyij_3,uzij_3 displacement at point i due to unit force in j
% in direction 3 (vertical one)
%Xi,Yi,Zi,Xj,Yj,Zj global coordinates of i and j point

%To account that Vaziri and Sun y axis is opposite to Midlin and my RF
Yi=-Yi;
Yj=-Yj;
patchB = -patchB;
patchT = -patchT;

%Soil properties
ni=nus;
Es=Gs.*2.*(1+ni);
beta=8.*pi.*Es.*(1-ni)./(1+ni);

%Coordinates of the pressure distribution location

% u1 = Xj-h_el_foot/2;
% u2 = Xj+h_el_foot/2;
% v1 = Yj-bfoot/2;
% v2 = Yj+bfoot/2;
u1 = patchL;
u2 = patchR;
v1 = patchT;
v2 = patchB;
w  = Zi;

%Pressure
p=1./(u2-u1)./(v2-v1);
q=1./(u2-u1)./(v2-v1);

%Coordinates of the poit at which the displacement due to the pressure is
%computed
x=Xi;
y=Yi;
z=Zi;



%-------------------------------------------------------------------
%Displacements due to vertical pressure p on a horizontal rectangle
%Eq(12) of Vaziri et al.
[ fdz_p_u1_v1 ] = f_int_mind_dz_p(u1,v1,w,x,y,z,ni);
[ fdz_p_u1_v2 ] = f_int_mind_dz_p(u1,v2,w,x,y,z,ni);
[ fdz_p_u2_v1 ] = f_int_mind_dz_p(u2,v1,w,x,y,z,ni);
[ fdz_p_u2_v2 ] = f_int_mind_dz_p(u2,v2,w,x,y,z,ni);
deltaz_p=p./beta.*(  (fdz_p_u2_v2-fdz_p_u1_v2)-(fdz_p_u2_v1-fdz_p_u1_v1));

%Eq(11) of Vaziri et al.
[ fdx_p_u1_v1 ] = f_int_mind_dx_p(u1,v1,w,x,y,z,ni);
[ fdx_p_u1_v2 ] = f_int_mind_dx_p(u1,v2,w,x,y,z,ni);
[ fdx_p_u2_v1 ] = f_int_mind_dx_p(u2,v1,w,x,y,z,ni);
[ fdx_p_u2_v2 ] = f_int_mind_dx_p(u2,v2,w,x,y,z,ni);
deltax_p=p./beta.*(  (fdx_p_u2_v2-fdx_p_u1_v2)-(fdx_p_u2_v1-fdx_p_u1_v1));

%Eq(11) of Vaziri et al.
[ fdy_p_u1_v1 ] = f_int_mind_dx_p(v1,u1,w,y,x,z,ni);
[ fdy_p_u1_v2 ] = f_int_mind_dx_p(v1,u2,w,y,x,z,ni);
[ fdy_p_u2_v1 ] = f_int_mind_dx_p(v2,u1,w,y,x,z,ni);
[ fdy_p_u2_v2 ] = f_int_mind_dx_p(v2,u2,w,y,x,z,ni);
deltay_p=p./beta.*(  (fdy_p_u2_v2-fdy_p_u1_v2)-(fdy_p_u2_v1-fdy_p_u1_v1));


%-------------------------------------------------------------------


uxij_3=deltax_p;
% uyij_3=-1*0.*y;
uyij_3=-deltay_p;
uzij_3=deltaz_p;





end


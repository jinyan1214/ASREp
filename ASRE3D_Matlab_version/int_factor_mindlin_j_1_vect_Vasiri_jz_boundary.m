function [ uxij_1,uyij_1,uzij_1 ] = int_factor_mindlin_j_1_vect_Vasiri_jz_boundary( Xi,Yi,Zi,Xj,Yj,Zj,Gs,nus,patchL, patchR, patchB, patchT) 
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

u1 = patchL*ones(size(Xj));
u2 = patchR*ones(size(Xj));
v1 = patchT*ones(size(Xj));
v2 = patchB*ones(size(Xj));


w  = Zi;

%Pressure
p=1./(u2-u1)./(v2-v1);
q=1./(u2-u1)./(v2-v1);

%Coordinates of the poit at which the displacement due to the pressure is
%computed
x=Xi;
y=Yi;
z=Zi;



%Displacements due to horizontal uniform shear stress q in direction x  on
%a horizontal rectangel 
%Eq(6) of Vaziri et al.
[ fdz_q_u1_v1 ] = f_int_mind_dz_q(u1,v1,w,x,y,z,ni);
[ fdz_q_u1_v2 ] = f_int_mind_dz_q(u1,v2,w,x,y,z,ni);
[ fdz_q_u2_v1 ] = f_int_mind_dz_q(u2,v1,w,x,y,z,ni);
[ fdz_q_u2_v2 ] = f_int_mind_dz_q(u2,v2,w,x,y,z,ni);
deltaz_q=q./beta.*(  (fdz_q_u2_v2-fdz_q_u1_v2)-(fdz_q_u2_v1-fdz_q_u1_v1));

%See scratch of 20210611
[ fdy_q_u1_v1 ] = f_int_mind_dy_q(u1,v1,w,x,y,z,ni);
[ fdy_q_u1_v2 ] = f_int_mind_dy_q(u1,v2,w,x,y,z,ni);
[ fdy_q_u2_v1 ] = f_int_mind_dy_q(u2,v1,w,x,y,z,ni);
[ fdy_q_u2_v2 ] = f_int_mind_dy_q(u2,v2,w,x,y,z,ni);
deltay_q=q./beta.*(  (fdy_q_u2_v2-fdy_q_u1_v2)-(fdy_q_u2_v1-fdy_q_u1_v1));

% %Eq(5) of Vaziri et al.
[ fdx_q_u1_v1 ] = f_int_mind_dx_q(u1,v1,w,x,y,z,ni);
[ fdx_q_u1_v2 ] = f_int_mind_dx_q(u1,v2,w,x,y,z,ni);
[ fdx_q_u2_v1 ] = f_int_mind_dx_q(u2,v1,w,x,y,z,ni);
[ fdx_q_u2_v2 ] = f_int_mind_dx_q(u2,v2,w,x,y,z,ni);
deltax_q=q./beta.*(  (fdx_q_u2_v2-fdx_q_u1_v2)-(fdx_q_u2_v1-fdx_q_u1_v1));


uxij_1=deltax_q;
uyij_1=-deltay_q;
uzij_1=deltaz_q;





end


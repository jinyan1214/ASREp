function [ term ] = f_int_mind_dx_p(u,v,w,x,y,z,ni )

%Eq(11) of Vaziri et al.

script_int_mind

term=-(3-4.*ni).*Z1.*log(Y+R2)...
    -Z1.*log(Y+R1)-(2.*w.*z.*Y.*Z2)./(R2.*(X.^2+Z2.^2))...
    -4.*(1-ni).*(1-2.*ni).*(Y.*log(Z2+R2)...
    +Z2.*log(Y+R2)+2.*X.*Tx2);
end


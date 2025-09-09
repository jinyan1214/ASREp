function [Talphai] = f_int_mind_Talphai(alpha,i,X,Y,R1,R2,Z1,Z2)

if i==1
    Zi=Z1;
    Ri=R1;
end

if i==2
    Zi=Z2;
    Ri=R2;
end

Talphai=atan((X+Y+Zi-alpha+Ri)./(alpha));



end


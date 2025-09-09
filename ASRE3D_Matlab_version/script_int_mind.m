X=x-u;
Y=y-v;
Z1=z-w;
Z2=z+w;
R1=(X.^2+Y.^2+Z1.^2).^0.5;
R2=(X.^2+Y.^2+Z2.^2).^0.5;

% Talphai=f_int_mind_Talphai(alpha,i,X,Y,R1,R2,Z1,Z2);
Tx1=f_int_mind_Talphai(X,1,X,Y,R1,R2,Z1,Z2);
Tx2=f_int_mind_Talphai(X,2,X,Y,R1,R2,Z1,Z2);
Ty1=f_int_mind_Talphai(Y,1,X,Y,R1,R2,Z1,Z2);
Ty2=f_int_mind_Talphai(Y,2,X,Y,R1,R2,Z1,Z2);
Tz1=f_int_mind_Talphai(Z1,1,X,Y,R1,R2,Z1,Z2);
Tz2=f_int_mind_Talphai(Z2,2,X,Y,R1,R2,Z1,Z2);

% Sun
Txy=f_int_mind_TAB(X,Y);
Tyx=f_int_mind_TAB(Y,X);

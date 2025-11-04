function A = A_mat(Y01,Y00,YG1,YG0,as,bs,at,bt,gs,gt,dt,ds,j)

global m I1 I3 ro
I=eye(3);

K01=K_mat(Y01(7:9),j);
KG1=K_mat(YG1(7:9),j);

K00=K_mat(Y00(7:9),j);
KG0=K_mat(YG0(7:9),j);
K=(1-bt)*((1-bs)*KG1+bs*KG0)+bt*((1-bs)*K01+bs*K00);

M=zeros(12,12);
M(4:6,7:9)=I;
M(7:9,4:6)=(ro*I1).*I;
M(9,6) = ro*I3;
M(10:12,1:3)=m*I;

A = ((1-at)*(1-as)/(gt*dt)).*M + ((1-bt)*((1-bs))/(gs*ds)).*K + ...
    ((1-bs)*(1-bt)).*jacobian_Fy(YG1,j);
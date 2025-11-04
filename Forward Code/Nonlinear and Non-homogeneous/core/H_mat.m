function H =H_mat(Ydot01,Ydot00,Yprime01,Yprime00,Y01,Y00,YG0,YG1,as,at,bs,bt,gs,gt,ds,dt,j,d)
global m I1 I3 ro

K01=K_mat(Y01(7:9),j);
KG1=K_mat(YG1(7:9),j);

I=eye(3);

K00=K_mat(Y00(7:9),j);
KG0=K_mat(YG0(7:9),j);
K=(1-bt)*((1-bs)*KG1+bs*KG0)+bt*((1-bs)*K01+bs*K00);

% M matrix

M=zeros(12,12);
M(4:6,7:9)=I;
M(7:9,4:6)=(ro*I1).*I;
M(9,6) = ro*I3;
M(10:12,1:3)=m*I;

H2 = ((1-at)*(1-as)*(1-gt)/gt - at*(1-as)).*(M*Ydot01) + ...
    ((1-at)*(as)*(1-gt)/gt - at*(as)).*(M*Ydot00) - ...
    (bt*(1-bs)).*(K*Yprime01) - (bt*bs).*(K*Yprime00)+...
    ((1-at)*(1-as)/(gt*dt)).*(M*Y01) + ...
    ((1-at)*(as)/(gt*dt)).*(M*Y00) - ...
    (bt*(1-bs)).*F_vec(Y01,j,d) - ...
    (bt*bs).*F_vec(Y00,j,d);

P = (bs-(1-bs)*(1-gs)/gs)*(1-bt);

YdotG = (YG0-Y00)/(gt*dt) - ((1-gt)/gt).*Ydot00;

H = H2 + P.*(F_vec(YG0,j,d)+M*YdotG) - ...
    ((1-bt)*bs).*F_vec(YG0,j,d)+...
    ((1-bt)*bs).*jacobian_Fy(YG0,j)*YG0...
    - ((1-bt)*(1-bs)).*F_vec(YG1,j,d)+...
    ((1-bs)*(1-bt)).*jacobian_Fy(YG1,j)*YG1;


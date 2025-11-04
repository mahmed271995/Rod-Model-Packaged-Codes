function  jac= jacobian_FX(X,j)
global E G m I1 I2 I3 ro
jac=zeros(12,12);

ds = 0.01;
I=eye(3);

C10=0;
C11=E*I1;
C12=0;
C13=0;
C14=0;
C15=0;
C20=0;
C21=E*I2;
C22=0;
C23=600;
C24=0;
C25=0;
C30=0;
C31=G*I3;
C32=0;
C33=0;
C34=0;
C35=0;

M=zeros(12,12);
M(4:6,7:9)=I;
M(7:9,4:6)=(ro*I1).*I;
M(9,6) = ro*I3;
M(10:12,1:3)=m*I;
% moment of inertias
Ixx=ro*I1;
Iyy=ro*I2;
Izz=ro*I3;
Im=[Ixx 0 0; 0 Iyy 0; 0 0 Izz];
% initial curvature:
k0=[0; 0; 0]; %(pi/32)

one_skew=-skew_conv([0;0;1]);   %[0 1 0;-1 0 0;0 0 0];
w=X(4:6,1);
k=X(7:9,1);

v_skew=skew_conv(X(1:3,1));
w_skew=skew_conv(X(4:6,1));
k_skew=skew_conv(X(7:9,1));
f_skew=skew_conv(X(10:12,1));

% wI_skew=[0 w(3,1)*(Izz-Iyy) w(2,1)*(Izz-Iyy);
%         -1*w(3,1)*(Izz-Ixx) 0 -1*w(1,1)*(Izz-Ixx);
%         w(2,1)*(Iyy-Ixx) w(1,1)*(Ixx-Iyy) 0 ];

%kCA_skew=[0 -1*k(3,1)*(C-A) -1*k(2,1)*(C-A); 
%          k(3,1)*(C-A) 0 k(1,1)*(C-A); 0 0 0];
%          kCA_skew=jackcl_mat(jackcl,X,ks1,ks2,ks3);

fs = 2*(ds*(j-1))^3-4*(ds*(j-1))^2+3*(ds*(j-1))+1;
fs_prime = 6*(ds*(j-1))^2-8*(ds*(j-1))+3;

k0=[0; 0; 0];

wI_skew=w_skew*Im - skew_conv(Im*w);

psi = [fs*(C10+C11*(X(7)-k0(1,1))+C12*((X(7)-k0(1,1))^2)+C13*((X(7)-k0(1,1))^3)+C14*((X(7)-k0(1,1))^4)+C15*((X(7)-k0(1,1))^5));
       fs*(C20+C21*(X(8)-k0(2,1))+C22*((X(8)-k0(2,1))^2)+C23*((X(8)-k0(2,1))^3)+C24*((X(8)-k0(2,1))^4)+C25*((X(8)-k0(2,1))^5));
       fs*(C30+C31*(X(9)-k0(3,1))+C32*((X(9)-k0(3,1))^2)+C33*((X(9)-k0(3,1))^3)+C34*((X(9)-k0(3,1))^4)+C35*((X(9)-k0(3,1))^5))];

psi_skew = skew_conv(psi);

dpsi_dk = [fs*(C11+2*C12*(X(7)-k0(1,1))+3*C13*((X(7)-k0(1,1))^2)+4*C14*((X(7)-k0(1,1))^3)+5*C15*((X(7)-k0(1,1))^4)) 0 0;
           0 fs*(C21+2*C22*(X(8)-k0(2,1))+3*C23*((X(8)-k0(2,1))^2)+4*C24*((X(8)-k0(2,1))^3)+5*C25*((X(8)-k0(2,1))^4)) 0;
           0 0 fs*(C31+2*C32*(X(9)-k0(3,1))+3*C33*((X(9)-k0(3,1))^2)+4*C34*((X(9)-k0(3,1))^3)+5*C35*((X(9)-k0(3,1))^4))];

Jpsi_s = -1*[fs_prime*(C11+2*C12*(X(7)-k0(1,1))+3*C13*((X(7)-k0(1,1))^2)+4*C14*((X(7)-k0(1,1))^3)+5*C15*((X(7)-k0(1,1))^4)) 0 0;
             0 fs_prime*(C21+2*C22*(X(8)-k0(2,1))+3*C23*((X(8)-k0(2,1))^2)+4*C24*((X(8)-k0(2,1))^3)+5*C25*((X(8)-k0(2,1))^4)) 0;
             0 0 fs_prime*(C31+2*C32*(X(9)-k0(3,1))+3*C33*((X(9)-k0(3,1))^2)+4*C34*((X(9)-k0(3,1))^3)+5*C35*((X(9)-k0(3,1))^4))];

kCA_skew = psi_skew - k_skew * dpsi_dk + Jpsi_s;

% kCA_skew=[0 k(3,1)*P1-R Q-k(2,1)*P1;
%           R-k(3,1)*Q1 0 k(1,1)*Q1-P;
%           k(2,1)*R1-Q P-k(1,1)*R1 0]';
      
jac=[-k_skew one_skew v_skew zeros(3,3);
    zeros(3,3) -k_skew w_skew zeros(3,3);
    zeros(3,3) wI_skew kCA_skew one_skew;
    m*w_skew -m*v_skew f_skew -k_skew];
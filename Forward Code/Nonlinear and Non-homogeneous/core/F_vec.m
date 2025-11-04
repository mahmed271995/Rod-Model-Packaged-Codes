function f_mat = F_vec(Y,j,d)
global E G m I1 I2 I3 ro Lee radi

eL = 1;       % has to be defined in 'YG_new_direct' as well
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

% initial curvature:
k0=[0; 0; 0]; %(pi/32)



%2*(ds*(j-1))^3-4*(ds*(j-1))^2+3*(ds*(j-1))+1;
%6*(ds*(j-1))^2-8*(ds*(j-1))+3;

fs = 2*(ds*(j-1))^3-4*(ds*(j-1))^2+3*(ds*(j-1))+1;
fs_prime = 6*(ds*(j-1))^2-8*(ds*(j-1))+3;

B=[fs*(C10+C11*(Y(7)-k0(1,1))+C12*((Y(7)-k0(1,1))^2)+C13*((Y(7)-k0(1,1))^3)+C14*((Y(7)-k0(1,1))^4)+C15*((Y(7)-k0(1,1))^5));
   fs*(C20+C21*(Y(8)-k0(2,1))+C22*((Y(8)-k0(2,1))^2)+C23*((Y(8)-k0(2,1))^3)+C24*((Y(8)-k0(2,1))^4)+C25*((Y(8)-k0(2,1))^5));
   fs*(C30+C31*(Y(9)-k0(3,1))+C32*((Y(9)-k0(3,1))^2)+C33*((Y(9)-k0(3,1))^3)+C34*((Y(9)-k0(3,1))^4)+C35*((Y(9)-k0(3,1))^5))];


  c=[0 0 1]';
  wCa3=cross(Y(4:6,1),c);
  kCv=cross(Y(7:9,1),Y(1:3,1));
    
  o1=wCa3-kCv;
  
  kCw=cross(Y(7:9,1),Y(4:6,1));
  t2= -1* kCw;
  
  IT=[ro*I1 0 0;0 ro*I2 0;0 0 ro*I3];
  ITw=IT*Y(4:6,1);
  wCITw=cross(Y(4:6,1),ITw);
  fCa3=cross(Y(10:12,1),c);
  kCBk=cross(Y(7:9),B);
  Qrc=0;
  
  %del_q_del_s=zeros(3,1);
  del_q_del_s = [fs_prime*(C10+C11*(Y(7)-k0(1,1))+C12*((Y(7)-k0(1,1))^2)+C13*((Y(7)-k0(1,1))^3)+C14*((Y(7)-k0(1,1))^4)+C15*((Y(7)-k0(1,1))^5));
                 fs_prime*(C20+C21*(Y(8)-k0(2,1))+C22*((Y(8)-k0(2,1))^2)+C23*((Y(8)-k0(2,1))^3)+C24*((Y(8)-k0(2,1))^4)+C25*((Y(8)-k0(2,1))^5));
                 fs_prime*(C30+C31*(Y(9)-k0(3,1))+C32*((Y(9)-k0(3,1))^2)+C33*((Y(9)-k0(3,1))^3)+C34*((Y(9)-k0(3,1))^4)+C35*((Y(9)-k0(3,1))^5))];
%   der=0;
%   del_q_del_s(1:2,1)=((E*pi*(D^3)/16)*der).*Y(7:8,1); 
%   del_q_del_s(3,1)=(G*pi*(D^3)/8)*der*Y(9,1);
  
  t3=wCITw+fCa3-kCBk-Qrc-del_q_del_s;
  
  mwCv=m*cross(Y(4:6,1),Y(1:3,1));
  kCf=cross(Y(7:9,1),Y(10:12,1));
 
%% Drag      
   %F_fluid=[0; 0; 0];
D=2*radi;
rof=1000;
Cn=0.1;
Ct=0.01;
v=Y(1:3,1);
t=[0;0;1];
VcrosT=cross(v,t); 
F_fluid=zeros(3,1);
%F_fluid=-0.5*rof*D*(Cn*norm(VcrosT,2)*cross(t,VcrosT) + pi*Ct*(dot(v,t)*norm(dot(v,t),2)*t) ); 
  

%   if Lee>eL*0.9995  
    F_dis=0.0001;
%       F_dis1=0.00001; 
%   else
%       F_dis1=0;
%   end

  Frc=[0; 0; 0;];
  
  f4=mwCv-kCf-Frc-F_fluid; % only Frc and F_fluid are external forces
  
  f_mat=[o1;t2;t3;f4];
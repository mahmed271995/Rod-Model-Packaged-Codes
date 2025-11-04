function Yprime_current = Yprime(Ydot_current,Y_current,d)
global m I1 I3 ro
I=eye(3);
N = length(Y_current(1,:));
Yprim=zeros(12,N);


for i = 1:N
    M=zeros(12,12);
    M(4:6,7:9)=I;
    M(7:9,4:6) = ro*I1.*I;
    M(9,6) = ro*I3;
    M(10:12,1:3) = m*I;

    K= K_mat(Y_current(7:9,i),i);
    F_mat = F_vec(Y_current(:,i),i,d);
    RHS = -F_mat - M*Ydot_current(:,i);
    Yprim(:,i) = K\RHS;
end
Yprime_current=Yprim;
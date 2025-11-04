function dy0_dt = Ydot0(Y0,ds)

N = length(Y0(1,:));

%% solving M as 9x9 matrix and by equation  MY_dot + KY_dash _F(Y) = 0;
%% solve v1,v2,...,k1,k2,k3 and now assume f1_dot,f2_dot,f3_dot =0 for our case
%% but if Y0 is not zero then case will be different
dy0_dt = zeros(12,N);


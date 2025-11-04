function skew_mat = skew_conv(X)
skew_mat=zeros(3,3);
skew_mat(1,2)= -X(3,1);
skew_mat(1,3)=  X(2,1);
skew_mat(2,1)=  X(3,1);
skew_mat(2,3)= -X(1,1);
skew_mat(3,1)= -X(2,1);
skew_mat(3,2)=  X(1,1);




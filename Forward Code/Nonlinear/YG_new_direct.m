%% Prescribe load
function YG1 = YG_new_direct(co_basis,Ypar,dt,d,Y)

eL = 1;
ds = 0.01;

N = eL/ds +1;
leftkeep=zeros(6,1);
coefficient=zeros(6,1);

%% Right hand boundary conditions
%B7-B9: kappa values (i.e: kappa1, kappa2, kappa3) which translate to
%moment q based on constitutive law
%B10-B12: force values (i.e: f1, f2 ,f3)

    %Bending loading
    BV7 = 0;
    BV8 = 4*sin(2*pi*d*dt);
    BV9 = 0;
    BV10 = 0;
    BV11 = 0;
    BV12 = 0;

    %Shear Loading
    % BV7 = 0;
    % BV8 = 0;
    % BV9 = 0;
    % BV10 = 15000;
    % BV11 = 0;
    % BV12 = 0;

    leftkeep(1,1) = BV7 - Ypar(7,end);
    leftkeep(2,1) = BV8 - Ypar(8,end);
    leftkeep(3,1) = BV9 - Ypar(9,end);
    leftkeep(4,1) = BV10 - Ypar(10,end);
    leftkeep(5,1) = BV11 - Ypar(11,end);
    leftkeep(6,1) = BV12 - Ypar(12,end);

    x(1:6,1:6)=co_basis(7:12,1:6,end);


    temp=x;
    coefficient=temp\leftkeep;


    for i = 1:N
        temp=co_basis(1:12,1:6,i);
        YG1(:,i) = Ypar(:,i)+temp*coefficient ;
    end
    
end


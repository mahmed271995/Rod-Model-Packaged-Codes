%% Left hand boundary conditions
function YP1 = leftbound(YG_previous)
    v1=0;
    v2=0;
    v3=0;
    w1=0;
    w2=0;
    w3=0;
    k1=YG_previous(7,1);
    k2=YG_previous(8,1);
    k3=YG_previous(9,1);
    f1=YG_previous(10,1);
    f2=YG_previous(11,1);
    f3=YG_previous(12,1);
    YP1=[v1 v2 v3 w1 w2 w3 k1 k2 k3 f1 f2 f3]';
end
function Y_time_march = run_numeric(N, ds, dt, nT, as, at, bs, bt, gs, gt)

global E I1 I2 I3 G ro Ac radi m poiss

basis = eye(12);
k0 = [0; 0; 0];  % Initial curvature

Y0 = zeros(12, N);
for i = 1:N
    Y0(7:9, i) = k0;
end
YG = Y0;
YG_new = zeros(12, N);
Ydot_initial = Ydot0(Y0, ds);
Y_time_march = zeros(12, N, nT+1);
Y_time_march(:,:,1) = Y0;

for d = 1:nT
    Yprime_previous = Yprime(Ydot_initial, Y0, d);
    e = 12;
    if mod(d,50) == 0
        d
    end 
    while e > 1e-8
        YG_previous = YG;
        co_basis = zeros(12, 6, N);
        Ypar = zeros(12, N);

        YP1 = leftbound(YG_previous);
        for i = 1:6
            co_basis(:, i, 1) = basis(:, i+6);
        end
        Ypar(:, 1) = YP1(:, 1);

        H2toN = zeros(12, N-1);
        A2toN = zeros(12, 12, N-1);
        B2toN = zeros(12, 12, N-1);

        for i = 1:N-1
            H2toN(:, i) = H_mat(Ydot_initial(:, i+1), Ydot_initial(:, i), ...
                Yprime_previous(:, i+1), Yprime_previous(:, i), Y0(:, i+1), ...
                Y0(:, i), YG(:, i), YG(:, i+1), as, at, bs, bt, gs, gt, ds, dt, i+1, d);
            A2toN(:,:,i) = A_mat(Y0(:, i+1), Y0(:, i), YG(:, i+1), YG(:, i), as, bs, at, bt, gs, gt, dt, ds, i+1);
            B2toN(:,:,i) = B_mat(Y0(:, i+1), Y0(:, i), YG(:, i+1), YG(:, i), as, bs, at, bt, gs, gt, dt, ds, i+1);
        end

        for j = 1:N-1
            for i = 1:6
                co_basis(:, i, j+1) = A2toN(:,:,j) \ (-B2toN(:,:,j) * co_basis(:, i, j));
            end
        end

        for j = 1:N-1
            Ypar(:, j+1) = A2toN(:,:,j) \ (H2toN(:, j) - B2toN(:,:,j) * Ypar(:, j));
        end

        YG_prev = YG;
        YG = YG_new_direct(co_basis, Ypar, dt, d, YG_prev);

        e = norm(YG - YG_previous);
    end

    Y_time_march(:,:,d+1) = YG;
    Ydot_initial = Ydot(Y_time_march(:,:,d+1), Y_time_march(:,:,d), Ydot_initial, gt, dt);
    Y0 = YG;
end
end

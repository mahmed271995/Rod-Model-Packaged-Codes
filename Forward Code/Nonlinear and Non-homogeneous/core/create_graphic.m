function create_graphic(video_name, Y_time_march, ds)
    % Create a video from Y_time_march showing rod motion
    % INPUTS:
    %   video_name       - filename as string (without extension)
    %   Y_time_march     - (12 x N x nT) array of rod state over time
    %   ds               - spatial discretization step

    [~, N, nT] = size(Y_time_march);

    x1 = zeros(nT, N);
    x2 = zeros(nT, N);
    x3 = zeros(nT, N);

    for i = 1:nT
        K_temp = Y_time_march(7:9,:,i);
        R = zeros(3, N);
        L = eye(3);

        for j = 1:N-1
            theta_v = ds * ((K_temp(:, j) + K_temp(:, j+1)) / 2);

            theta_m = zeros(3);
            theta_m(1,2) = -theta_v(3);
            theta_m(1,3) =  theta_v(2);
            theta_m(2,1) =  theta_v(3);
            theta_m(2,3) = -theta_v(1);
            theta_m(3,1) = -theta_v(2);
            theta_m(3,2) =  theta_v(1);

            L = expm(-theta_m) * L;
            r = L(3,:)';
            R(:, j+1) = R(:, j) + r * ds;
        end

        x1(i,:) = R(1,:);
        x2(i,:) = R(2,:);
        x3(i,:) = R(3,:);
    end

    set(gcf,'Units','inches','Position',[0 1 16 7])
    set(0,'defaultTextInterpreter','latex')
    set(gcf,'color','w')
    v = VideoWriter([video_name, '.mp4'], 'MPEG-4');
    open(v);

    figure;
    for i = 1:5:nT

        subplot(2,2,1)
            plot(x3(i,:), x1(i,:), 'k', 'LineWidth', 2);
            set(gca, 'FontSize', 16);
            xlim([-0.5 1.5]);
            ylim([-1 1]);
            xlabel('$z$ (m)', 'FontSize', 22);
            ylabel('$x$ (m)', 'FontSize', 22);
            drawnow;
        subplot(2,2,2)
            plot(x3(i,:), x2(i,:), 'k', 'LineWidth', 2);
            set(gca, 'FontSize', 16);
            xlim([-0.5 1.5]);
            ylim([-1 1]);
            xlabel('$z$ (m)', 'FontSize', 22);
            ylabel('$y$ (m)', 'FontSize', 22);

            drawnow;
        subplot(2,2,3)
            plot(x1(i,:), x2(i,:), 'k', 'LineWidth', 2);
            set(gca, 'FontSize', 16);
            xlim([-1 1]);
            ylim([-1 1]);
            xlabel('$x$ (m)', 'FontSize', 22);
            ylabel('$y$ (m)', 'FontSize', 22);
            drawnow;
        subplot(2,2,4)
            plot3(x3(i,:),x1(i,:), x2(i,:), 'k', 'LineWidth', 2);
            set(gca, 'FontSize', 16);
            xlim([-0.5 1.5]);
            ylim([-1 1]);
            zlim([-1 1]);
            xlabel('$z$ (m)', 'FontSize', 22);
            ylabel('$x$ (m)', 'FontSize', 22);
            zlabel('$y$ (m)', 'FontSize', 22);
        frame = getframe(gcf);
        writeVideo(v, frame);
    end

    close(v);
    disp(['Video saved as ', video_name, '.mp4']);
end

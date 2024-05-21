function [repro_error, x2_repro] = backprojection(correspondences, P1, Image2, T, R, K)
    % This function calculates the mean error of the back projection
    % of the world coordinates P1 from image 1 in camera frame 2
    % and visualizes the correct feature coordinates as well as the back projected ones.
    N = size(P1, 2);
    x2 = correspondences(3:4, :); %2x29
    x2 = cat(1, x2, ones(1, size(x2, 2)));%3x29
    P1_2 = R*P1 + T; %3x29
    P1_2 = P1_2(:, :) ./ P1_2(3,:); %Scaling Z = 1
    x2_repro = K * P1_2; %3x29
    figure();
    hold all;
    imshow(Image2);
    repro_error = 0
    for i=1:N
        plot(x2(1,i), x2(2,i), 'bx');
        text(x2(1,i), x2(2,i), int2str(i),'color', 'blue');
        plot(x2_repro(1,i), x2_repro(2,i), 'ro');
        text(x2_repro(1,i), x2_repro(2,i), int2str(i),'color', 'red');
        repro_error = repro_error + norm(x2(:,i) - x2_repro(:,i), 2);
    end
    repro_error = repro_error / i;
end
function [T, R, lambda, P1, camC1, camC2] = reconstruction(T1, T2, R1, R2, correspondences, K)
    % This function estimates the depth information and thereby determines the 
    % correct Euclidean movement R and T. Additionally it returns the
    % world coordinates of the image points regarding image 1 and their depth information.
    
    %% Preparation from task 4.2
    % T_cell    cell array with T1 and T2 
    % R_cell    cell array with R1 and R2
    % d_cell    cell array for the depth information
    % x1        homogeneous calibrated coordinates
    % x2        homogeneous calibrated coordinates
    preparation
    
    %% R, T and lambda from task 4.3
    % T         reconstructed translation
    % R         reconstructed rotation
    % lambda    depth information
    R_T_lambda
    
    %% Calculation and visualization of the 3D points and the cameras
    P1 = lambda(:,1)'.*x1;
    camC1 = [-0.2, 0.2, 0.2, -0.2;
        0.2, 0.2, -0.2, -0.2;
        1, 1, 1, 1];
    Euk_Trafo = [R, T;
        zeros(1,3), 1];
    camC2 = Euk_Trafo \ [camC1;ones(1,4)];
    camC2 = camC2(1:end-1, :);
    camC1 = cat(2,camC1, camC1(:,1));
    camC2 = cat(2,camC2, camC2(:,1));
    figure();
    hold all;
    number = int2str(1:size(x1,2));
    scatter3(P1(1,:), P1(2,:), P1(3,:));
    text(P1(1,:), P1(2,:), P1(3,:), number);
    plot3(camC1(1,:), camC1(2,:), camC1(3,:),'b');
    plot3(camC2(1,:), camC2(2,:), camC2(3,:),'r');
    text(camC1(1,1), camC1(2,1), camC1(3,1), 'Cam1', 'Color','blue');
    text(camC2(1,1), camC2(2,1), camC2(3,1), 'Cam2', 'Color','red');
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
    camera_position = [43, -22, -87];    % [x, y, z] coordinates of the camera position
    camera_target = [0, -1, 0];          % [x, y, z] coordinates of the camera target
    campos(camera_position)
    camup(camera_target)
    camC1 = camC1(:, 1:end-1);
    camC2 = camC2(:, 1:end-1);
end
function [T, R, lambda, M1, M2] = reconstruction(T1, T2, R1, R2, correspondences, K)
    %% Preparation from task 4.2
    % T_cell    cell array with T1 and T2 
    % R_cell    cell array with R1 and R2
    % d_cell    cell array for the depth information
    % x1        homogeneous calibrated coordinates
    % x2        homogeneous calibrated coordinates
    preparation
    
    %% Reconstruction
    N = size(x1,2);
    amount_pos = 0;
    for i = 1:4
        T = T_cell{i};
        R = R_cell{i};
        M1 = zeros(3*N,N+1);
        M2 = zeros(3*N,N+1);
        for j=1:N
            M1(3*j-2:3*j,j) = hat(x2(:,j)) * R * x1(:,j);
            M1(3*j-2:3*j,N+1) = hat(x2(:,j)) * T;
            M2(3*j-2:3*j,j) = hat(x1(:,j)) * R' * x2(:,j);
            M2(3*j-2:3*j,N+1) = -hat(x1(:,j)) * R' * T;
        end
        [U1 S1 V1] = svd(M1);
        [U2 S2 V2] = svd(M2);
        d1 = V1(:,end);
        d2 = V2(:,end);
        d1 = d1(1:end-1)/d1(end);
        d2 = d2(1:end-1)/d2(end);
        d_cell{i} = [d1, d2];
        mask = find([d1, d2] > 0);
        if length(mask) > amount_pos
            amount_pos = length(mask(:));
            save_i = i
            lambda = [d1 d2];
            M_last1 = M1;
            M_last2 = M2;
        end
    end
    T = T_cell{save_i};
    R = R_cell{save_i};
end
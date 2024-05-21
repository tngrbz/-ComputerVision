function [T_cell, R_cell, d_cell, x1, x2] = reconstruction(T1, T2, R1, R2, correspondences, K)
    %% Preparation
    x1 = correspondences(1:2, :);
    x2 = correspondences(3:4, :);
    x1 = cat(1, x1, ones(1, size(x1, 2)));
    x2 = cat(1, x2, ones(1, size(x2, 2)));
    x1 = K\x1;
    x2 = K\x2;
    N = size(x1, 2);
    d_cell = cell(1,4);
    for i=1:4
        d_cell{i} = zeros(N,2);
    end
    T_cell = cell(1,4);
    R_cell = cell(1,4);
    T_cell{1} = T1;
    T_cell{2} = T2;
    T_cell{3} = T1;
    T_cell{4} = T2;
    R_cell{1} = R1;
    R_cell{2} = R1;
    R_cell{3} = R2;
    R_cell{4} = R2;
end
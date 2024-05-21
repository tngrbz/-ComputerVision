function [EF] = epa(correspondences, K)
    % Depending on whether a calibrating matrix 'K' is given,
    % this function calculates either the essential or the fundamental matrix
    % with the eight-point algorithm.
    
    %% First step of the eight-point algorithm from task 3.1
    % Known variables:
    % x1, x2        homogeneous (calibrated) coordinates       
    % A             matrix A for the eight-point algorithm
    % V             right-singular vectors
    epa_part1
    
    %% Estimation of the matrices
    g = V(:, end);
    G = reshape(g, [3, 3]);
    [Ug Sg Vg] = svd(G);
    
    if nargin == 1
        Sg(end, end) = 0;
        F = Ug * Sg * Vg';
        EF = F;
    else if nargin == 2
        Sg = eye(3);
        Sg(end, end) = 0;
        E = Ug * Sg * Vg';
        EF = E;
    end
end
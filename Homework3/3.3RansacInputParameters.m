function [epsilon, p, tolerance, x1_pixel, x2_pixel] = F_ransac(correspondences, varargin)
    % This function implements the RANSAC algorithm to determine 
    % robust corresponding image points
    defaultEps = 0.5;
    defaultP = 0.5;
    defaultTol = 0.01;
    
    p = inputParser;
    cond = @(x) isnumeric(x) && x>0 && x<1
    addRequired(p, 'correspondences', @(x) size(x, 1) > 0);
    addParameter(p, 'epsilon', defaultEps, cond);
    addParameter(p, 'p', defaultP, cond);
    addParameter(p, 'tolerance', defaultTol, @(x) isnumeric(x));
    parse(p, correspondences, varargin{:});
    epsilon = p.Results.epsilon;
    tolerance = p.Results.tolerance;
    p = p.Results.p;
    
    x1_pixel = correspondences(1:2, :);
    x2_pixel = correspondences(3:4, :);
    x1_pixel = cat(1, x1_pixel, ones(1, size(x1_pixel, 2)));
    x2_pixel = cat(1, x2_pixel, ones(1, size(x2_pixel, 2)));
end
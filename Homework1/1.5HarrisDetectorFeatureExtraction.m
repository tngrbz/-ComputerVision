function [H, corners, features] = harris_detector(input_image, varargin)
    % In this function you are going to implement a Harris detector that extracts features
    % from the input_image.
    
    %% Input parser from task 1.3
    % segment_length    size of the image segment
    % k                 weighting between corner- and edge-priority
    % tau               threshold value for detection of a corner
    % do_plot           image display variable
    input_parser

    %% Preparation for feature extraction from task 1.4
    % Ix, Iy            image gradient in x- and y-direction
    % w                 weighting vector
    % G11, G12, G22     entries of the Harris matrix
    image_preprocessing
    
    %% Feature extraction with the Harris measurement
    H = zeros(size(input_image));
    corners = zeros(size(input_image));
    for i = 1:size(H, 1)
        for j = 1:size(H, 2)
            A = G11(i, j);
            B = G12(i, j);
            C = G22(i, j);
            G = [A B; B C];
            R = det(G) - k * trace(G)^2;
            H(i, j) = R;
        end
    end
    cond = H > tau;
    corners(cond) = H(cond);
    % Eliminate corners close to the border
    segment_width = segment_length;
    border = ceil(segment_width/2);
    corners(1:border,:) = 0;
    corners(:,1:border) = 0;
    corners(end-border+1:end,:) = 0;
    corners(:,end-border+1:end) = 0;
    
    % Create matrix of corner coordinates
    [y,x] = find(corners);
    features = [x,y]';
end
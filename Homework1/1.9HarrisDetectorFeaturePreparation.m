function [corners, sorted_index] = harris_detector(input_image, varargin)
    % In this function you are going to implement a Harris detector that extracts features
    % from the input_image.
    
    %% Input parser from task 1.7
    % segment_length    size of the image segment
    % k                 weighting between corner- and edge-priority
    % tau               threshold value for detection of a corner
    % do_plot           image display variable
    % min_dist          minimal distance of two features in pixels
    % tile_size         size of the tiles
    % N                 maximal number of features per tile
    input_parser_new

    %% Preparation for feature extraction from task 1.4
    % Ix, Iy            image gradient in x- and y-direction
    % w                 weighting vector
    % G11, G12, G22     entries of the Harris matrix
    image_preprocessing
    
    %% Feature extraction with the Harris measurement from task 1.5
    % corners           matrix containing the value of the Harris measurement for each pixel         
    % features          detected features
    harris_measurement
    
    %% Feature preparation
    
    %Zero Border
    buffer = zeros(min_dist, size(corners, 2));
    corners = cat(1, corners, buffer);
    %size(corners)
    corners = cat(1, buffer, corners);
    %size(corners)
    buffer = zeros(size(corners, 1), min_dist);
    %size(buffer)
    corners = cat(2, corners, buffer);
    corners = cat(2, buffer, corners);
    
    % Sorting
    vector_corners = corners(:);
    non_zero_indices = find(vector_corners);
    vector_corners = vector_corners(non_zero_indices);
    [sorted_corner indices] = sort(vector_corners, 'descend');
    % [rows, cols] = ind2sub(size(corners), indices); % unravel index
    % sorted_index = [cols, rows];
    sorted_index = non_zero_indices(indices);
end
function [acc_array, features] = harris_detector(input_image, varargin)
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
    harris_measurement
    
    %% Feature preparation from task 1.9
    %corners            Harris measurement for each pixel respecting the minimal distance
    %sorted_index       Index list of features sorted descending by thier strength
    feature_preprocessing
    
    %% Accumulator array
    amount_tilesy = ceil(size(input_image, 1) / tile_size(1));
    amount_tilesx = ceil(size(input_image, 2) / tile_size(2));
    acc_array = zeros(amount_tilesy, amount_tilesx);
    features = zeros(2, N * size(acc_array(:), 1));
    if size(sorted_index, 1) < N * size(acc_array(:), 1)
        features = zeros(2, size(sorted_index, 1));
    end
end
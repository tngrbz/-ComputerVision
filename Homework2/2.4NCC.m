function [NCC_matrix, sorted_index] = point_correspondence(I1, I2, Ftp1, Ftp2, varargin)
    % In this function you are going to compare the extracted features of a stereo recording
    % with NCC to determine corresponding image points.
    
    %% Input parser from task 2.1
    % window_length         side length of quadratic window
    % min_corr              threshold for the correlation of two features
    % do_plot               image display variable
    % Im1, Im2              input images (double)
    input_parser
    
    %% Feature preparation from task 2.2
    % no_pts1, no_pts 2     number of features remaining in each image
    % Ftp1, Ftp2            preprocessed features
    feature_preprocessing
    
    %% Normalization from task 2.3
    % Mat_feat_1            normalized windows in image 1
    % Mat_feat_2            normalized windows in image 2
    window_normalization
    
    %% NCC calculations
    [y1, x1] = size(Mat_feat_1);
    [y2, x2] = size(Mat_feat_2);
    NCC_matrix = zeros(x2, x1);
    for i=1:x2
        W1 = reshape(Mat_feat_2(:, i), [window_length, window_length]);
        for j=1:x1
            W2 = reshape(Mat_feat_1(:, j), [window_length, window_length]);
            NCC_matrix(i,j) = 1/(window_length^2-1) * trace(W1'*W2);
        end
    end
    mask = NCC_matrix < min_corr;
    NCC_matrix(mask) = 0;
    mask = find(NCC_matrix(:));
    [~, idx] = sort(NCC_matrix(mask), 'descend');
    sorted_index = mask(idx);
end
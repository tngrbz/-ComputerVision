function cor = point_correspondence(I1, I2, Ftp1, Ftp2, varargin)
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
    
    %% NCC from task 2.4
    % NCC_matrix            matrix containing the correlation between the image points
    % sorted_index          sorted indices of NCC_matrix entries in decreasing order of intensity
    ncc_calculation
    
    %% Correspondeces
    [y, x] = ind2sub(size(NCC_matrix), sorted_index);
    cor = [0 0 0 0]';
    for i=1:size(sorted_index, 1)
        x_ = x(i);
        y_ = y(i);
        if NCC_matrix(y_, x_) ~= 0   
            NCC_matrix(:, x_) = zeros(size(NCC_matrix, 1), 1);
            coor1 = Ftp1(:, x_);
            coor2 = Ftp2(:, y_);
            cor = cat(2, cor, [coor1; coor2]);
        end
    end
    cor = cor(:, 2:end);
end
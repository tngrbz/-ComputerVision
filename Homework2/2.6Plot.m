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
    
    %% Correspondeces from task 2.5
    % cor                   matrix containing all corresponding image points
    correspondence
    
    %% Visualize the correspoinding image point pairs
    if do_plot
        im1_xy = cor(1:2, :);
        im2_xy = cor(3:4, :);
        
        figure(1);
        imshow(I1);
        hold on;
        im = imshow(I2);
        set(im, 'AlphaData',  0.5);
        plot(im1_xy(1,:), im1_xy(2,:), 'o');
        plot(im2_xy(1,:), im2_xy(2,:), 'o');
        for i=1:size(cor, 2)
            plot([im1_xy(1,i), im2_xy(1,i)], [im1_xy(2,i), im2_xy(2,i)]);
        end
end
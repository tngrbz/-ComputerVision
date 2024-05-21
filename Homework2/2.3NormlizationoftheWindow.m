function [Mat_feat_1, Mat_feat_2] = point_correspondence(I1, I2, Ftp1, Ftp2, varargin)
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
    
    %% Normalization
    Mat_feat_1 = zeros(window_length^2, 1);
    for i=1:size(Ftp1, 2)
        x = Ftp1(1, i);
        y = Ftp1(2, i);
        exp = (window_length-1)/2;
        W = Im1(y-exp:y+exp, x-exp:x+exp);
        N = size(W(:),1);
        avgW = 1/N*ones(size(W))*W*ones(size(W));
        std_dvtn = sqrt(norm(W-avgW, "fro")^2/(N-1));
        normW = (W-avgW)/std_dvtn;
        Mat_feat_1 = cat(2, Mat_feat_1, normW(:));
    end
    Mat_feat_1 = Mat_feat_1(:,2:end);
    
    %Image 2
    Mat_feat_2 = zeros(window_length^2, 1);
    for i=1:size(Ftp2, 2)
        x = Ftp2(1, i);
        y = Ftp2(2, i);
        exp = (window_length-1)/2;
        W = Im2(y-exp:y+exp, x-exp:x+exp);
        N = size(W(:),1);
        avgW = 1/N*ones(size(W))*W*ones(size(W));
        std_dvtn = sqrt(norm(W-avgW, "fro")^2/(N-1));
        normW = (W-avgW)/std_dvtn;
        Mat_feat_2 = cat(2, Mat_feat_2, normW(:));
    end
    Mat_feat_2 = Mat_feat_2(:,2:end);
end
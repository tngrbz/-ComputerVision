function [no_pts1, no_pts2, Ftp1, Ftp2] = point_correspondence(I1, I2, Ftp1, Ftp2, varargin)
    % In this function you are going to compare the extracted features of a stereo recording
    % with NCC to determine corresponding image points.
    
    %% Input parser from task 2.1
    % window_length     side length of quadratic window
    % min_corr          threshold for the correlation of two features
    % do_plot           image display variable
    % Im1, Im2          input images (double)
    input_parser
    
    %% Feature preparation
    % Image 1
    Im1_size = size(Im1);
    size(Ftp1)
    for i=1:size(Ftp1, 2)
        if Ftp1(1,i) <= (window_length-1)/2
            Ftp1(:,i) = [0;0];
        elseif Im1_size(2) - Ftp1(1,i) <= (window_length-1)/2
            Ftp1(:,i) = [0;0];
        elseif Ftp1(2,i) <= (window_length-1)/2
            Ftp1(:,i) = [0;0];
        elseif Im1_size(1) - Ftp1(2,i) <= (window_length-1)/2
            Ftp1(:,i) = [0;0];
        end
    end
    idx = find(Ftp1(1,:));
    Ftp1 = Ftp1(:,idx);
    no_pts1 = size(Ftp1,2);
    
    % Image 2
    Im2_size = size(Im2);
    size(Ftp2)
    for i=1:size(Ftp2, 2)
        if Ftp2(1,i) <= (window_length-1)/2
            Ftp2(:,i) = [0;0];
        elseif Im2_size(2) - Ftp2(1,i) <= (window_length-1)/2
            Ftp2(:,i) = [0;0];
        elseif Ftp2(2,i) <= (window_length-1)/2
            Ftp2(:,i) = [0;0];
        elseif Im2_size(1) - Ftp2(2,i) <= (window_length-1)/2
            Ftp2(:,i) = [0;0];
        end
    end
    idx = find(Ftp2(1,:));
    Ftp2 = Ftp2(:,idx);
    no_pts2 = size(Ftp2,2);
end
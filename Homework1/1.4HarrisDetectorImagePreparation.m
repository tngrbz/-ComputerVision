function [Ix, Iy, w, G11, G22, G12] = harris_detector(input_image, varargin)
    % In this function you are going to implement a Harris detector that extracts features
    % from the input_image.
    
    %% Input parser from task 1.3
    % segment_length    size of the image segment
    % k                 weighting between corner- and edge-priority
    % tau               threshold value for detection of a corner
    % do_plot           image display variable
    input_parser
    
    %% Preparation for feature extraction
    % Check if it is a grayscale image
    if size(input_image, 3) ~= 1
        error("Image format has to be NxMx1")
    end
    % Approximation of the image gradient
    [Ix, Iy] = sobel_xy(input_image);
    % Weighting
    segment_width = segment_length;

    % Calculate standard deviation of Gaussian distribution
    sigma = segment_width/6;
    
    % Create Gaussian weighting vector
    x = linspace(-segment_width/2, segment_width/2, segment_width);
    w = exp(-x.^2/(2*sigma^2));
    
    % Normalize weighting vector
    w = w/sum(w);
    W = w' .* w;
    % Harris Matrix G
    G11 = conv2(Ix .^2, W, 'same');
    G12 = conv2(Ix .* Iy, W, 'same');
    G22 = conv2(Iy .^2, W, 'same');
            
end
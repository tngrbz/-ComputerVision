function [segment_length, k, tau, do_plot] = harris_detector(input_image, varargin)
    % In this function you are going to implement a Harris detector that extracts features
    % from the input_image.

    %% Input parser
    segment_length = 15;
    k = 0.05;
    tau = 1000000;
    do_plot = false;
    
    p = inputParser;
    addRequired(p, 'input_image');
    addParameter(p, 'segment_length', segment_length, @(x) isnumeric(x) && (mod(x, 2) ~= 0) && (x>1));
    addParameter(p, 'k', k, @(x) isnumeric(x) && (x>=0) && (x<=1));
    addParameter(p, 'tau', tau, @(x) (x > 0));
    addParameter(p, 'do_plot', do_plot, @(x) islogical(x));
    parse(p,input_image,varargin{:});
    segment_length = p.Results.segment_length;
    k = p.Results.k;
    tau = p.Results.tau;
    do_plot = p.Results.do_plot;
    input_image = p.Results.input_image;
end
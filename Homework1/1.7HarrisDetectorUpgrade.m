function [min_dist, tile_size, N] = harris_detector(input_image, varargin)
    % In this function you are going to implement a Harris detector that extracts features
    % from the input_image.
    
    segment_length = 15;
    k = 0.05;
    tau = 1000000;
    do_plot = false;
    min_dist = 20;
    tile_size = [200, 200];
    N = 5;
    %% Input parser
    p = inputParser;
    addRequired(p, 'input_image');
    addParameter(p, 'segment_length', segment_length, @(x) isnumeric(x) && (mod(x, 2) ~= 0) && (x>1));
    addParameter(p, 'k', k, @(x) isnumeric(x) && (x>=0) && (x<=1));
    addParameter(p, 'tau', tau, @(x) (x > 0));
    addParameter(p, 'do_plot', do_plot, @(x) islogical(x));
    addParameter(p, 'min_dist', min_dist, @(x) (x>=1) && isnumeric(x));
    addParameter(p, 'tile_size', tile_size, @(x) isnumeric(x) && size(x,2) <= 2);
    addParameter(p, 'N', N, @(x) (x>=5) && isnumeric(x));
    
    parse(p,input_image,varargin{:});
    segment_length = p.Results.segment_length;
    k = p.Results.k;
    tau = p.Results.tau;
    do_plot = p.Results.do_plot;
    min_dist = p.Results.min_dist;
    tile_size = p.Results.tile_size;
    if size(tile_size, 2) == 1
        tile_size = [tile_size, tile_size];
    end
    N = p.Results.N;
    input_image = p.Results.input_image;
end
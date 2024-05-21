function [window_length, min_corr, do_plot, Im1, Im2] = point_correspondence(I1, I2, Ftp1, Ftp2, varargin)
    % In this function you are going to compare the extracted features of a stereo recording
    % with NCC to determine corresponding image points.
    
    %% Input parser
    Im1 = double(I1);
    Im2 = double(I2);
    defaultWindowLength = 25;
    defaultMinCorr = 0.95;
    defaultPlot = false;
    
    p = inputParser;
    addRequired(p, 'I1');
    addRequired(p, 'I2');
    addRequired(p, 'Ftp1');
    addRequired(p, 'Ftp2');
    addParameter(p, 'window_length', defaultWindowLength,@(x) isnumeric(x) && mod(x, 2) == 1 && x > 1);
    addParameter(p, 'min_corr', defaultMinCorr, @(x) isnumeric(x) && x > 0 && x < 1);
    addParameter(p, 'do_plot', defaultPlot, @(x) islogical(x));
    
    parse(p, I1, I2, Ftp1, Ftp2, varargin{:});
    window_length = p.Results.window_length;
    min_corr = p.Results.min_corr;
    do_plot = p.Results.do_plot;

end
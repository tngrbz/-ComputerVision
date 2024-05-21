function features = harris_detector(input_image, varargin)
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
    % sorted_index      sorted indices of features in decreasing order of feature strength
    feature_preprocessing
    
    %% Accumulator array from task 1.10
    % acc_array         accumulator array which counts the features per tile
    % features          empty array for storing the final features
    accumulator_array
    
    %% Feature detection with minimal distance and maximal number of features per tile
    Cake = cake(min_dist);
    Cake(min_dist + 1, min_dist + 1) = 1;
    for i = 1:size(sorted_index, 1)
        index = sorted_index(i);
        [rows, cols] = ind2sub(size(corners), index);
        if corners(rows, cols) ~= 0
            buffer = corners(rows - min_dist: rows + min_dist, cols - min_dist: cols + min_dist) .* Cake;
            size(buffer);
            min_dist;
            corners(rows - min_dist: rows + min_dist, cols - min_dist: cols + min_dist) = buffer;
        end
    end
    features = [0;0];
    tiles = cell(size(acc_array));
    corners = corners(min_dist+1:end-min_dist, min_dist+1:end-min_dist);
    for y = 1:tile_size(1):size(corners, 1)
        if y-1+tile_size(1) > size(corners, 1)
            upper_rangey = size(corners, 1);
        else
            upper_rangey = y-1+tile_size(1);
        end
        for x = 1:tile_size(2):size(corners, 2)
            if x-1+tile_size(2) > size(corners, 2)
            upper_rangex = size(corners, 2);
            else
            upper_rangex = x-1+tile_size(2);
            end
            tiles{floor(y/tile_size(1)) + 1, floor(x/tile_size(2)) + 1} = corners(y:upper_rangey, x:upper_rangex);
            %{
            vector_tile = corners(y:upper_rangey, x:upper_rangex);
            vector_tile = vector_tile(:);
            non_zero_indices = find(vector_tile);
            vector_tile = vector_tile(non_zero_indices);
            [~, sorted] = sort(vector_tile, 'descend');
            if size(sorted, 1) > N
                sorted = sorted(1:N, 1);
            end
            sorted = non_zero_indices(sorted);
            [rows, cols] = ind2sub(size(tiles{i,j}), sorted);
            %}
        end
    end
    
    for j=1:size(acc_array, 2)
        for i=1:size(acc_array, 1)
            vector_tile = tiles{i,j};
            vector_tile = vector_tile(:);
            non_zero_indices = find(vector_tile);
            vector_tile = vector_tile(non_zero_indices);
            [~, sorted] = sort(vector_tile, 'descend');
            if size(sorted, 1) > N
                sorted = sorted(1:N, 1);
            end
            sorted = non_zero_indices(sorted);
            [rows, cols] = ind2sub(size(tiles{i,j}), sorted);
            rows_shift = ((i-1) * tile_size(1));
            cols_shift = ((j-1) * tile_size(2));
            rows = rows_shift + rows;
            cols = cols_shift + cols;
            features = cat(2, features, [cols'; rows']);
        end
    end
    features = features(:, 2:end);
    %features = features - [min_dist;min_dist];
    size(features)
    
    % Plot Routine
    plotting
    
end
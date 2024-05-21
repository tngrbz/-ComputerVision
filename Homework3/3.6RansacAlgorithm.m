function [correspondences_robust, largest_set_F] = F_ransac(correspondences, varargin)
    % This function implements the RANSAC algorithm to determine 
    % robust corresponding image points
       
    %% Input parser
    % Known variables:
    % epsilon       estimated probability
    % p             desired probability
    % tolerance     tolerance to belong to the consensus-set
    % x1_pixel      homogeneous pixel coordinates
    % x2_pixel      homogeneous pixel coordinates
    input_parser
        
    %% RANSAC algorithm preparation
    % Pre-initialized variables:
    % k                     number of necessary points
    % s                     iteration number
    % largest_set_size      size of the so far biggest consensus-set
    % largest_set_dist      Sampson distance of the so far biggest consensus-set
    % largest_set_F         fundamental matrix of the so far biggest consensus-set
    ransac_preparation
    
    %% RANSAC algorithm
    x1_pre = correspondences(1:2, :);
    x2_pre = correspondences(3:4, :);
    x1 = cat(1, x1_pre, ones(1, size(x1_pre, 2)));
    x2 = cat(1, x2_pre, ones(1, size(x2_pre, 2)));
    sz_cor = size(correspondences, 2);
    correspondences_robust = [];
    for i=1:s
        rand = randi(sz_cor, [1 k]);
        F = epa(correspondences(:, rand));
        sd = sampson_dist(F, x1, x2);
        mask = sd<tolerance;
        cons = sd(mask);
        amount_cons = length(cons);
        dist_cons = sum(cons, 2);
        rob_x1 = x1_pre(:, mask);
        rob_x2 = x2_pre(:, mask);
        %if dist_cons > largest_set_distance
        %    largest_set_distance = dist_cons;
        %end
        if amount_cons > largest_set_size
            largest_set_size = amount_cons;
            largest_set_F = F;
            largest_set_distance = dist_cons;
            correspondences_robust = [rob_x1; rob_x2];
        elseif amount_cons == largest_set_size
            if dist_cons < largest_set_distance
                largest_set_size = amount_cons;
                largest_set_F = F;
                largest_set_distance = dist_cons;
                correspondences_robust = [rob_x1; rob_x2];
            end
        end
    end

end
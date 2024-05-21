function Cake = cake(min_dist)
    % The cake function creates a "cake matrix" that contains a circular set-up of zeros
    % and fills the rest of the matrix with ones. 
    % This function can be used to eliminate all potential features around a stronger feature
    % that don't meet the minimal distance to this respective feature.
    Cake = zeros(min_dist * 2 + 1, min_dist * 2 + 1);
    center = Cake(min_dist+1, min_dist+1);
    for i = 1:size(Cake, 1)
        for j = 1:size(Cake, 2)
            if sqrt((i - min_dist-1) ^2 + (j - min_dist-1)^2) > min_dist
                Cake(i,j) = 1;
            end
        end
    end
    Cake = logical(Cake);
end
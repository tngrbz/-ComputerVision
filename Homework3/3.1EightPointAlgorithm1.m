function [x1, x2, A, V] = epa(correspondences, K)
    % Depending on whether a calibrating matrix 'K' is given,
    % this function calculates either the essential or the fundamental matrix
    % with the eight-point algorithm.
    if nargin == 1
        K = [];
        disp('K is not given!')
    elseif nargin ==2
        disp('K is given!')
    end
    
    x1 = correspondences(1:2, :);
    x2 = correspondences(3:4, :);
    x1 = cat(1, x1, ones(1, size(x1, 2)));
    x2 = cat(1, x2, ones(1, size(x2, 2)));
    if ~isempty(K)
        x1 = K\x1;
        x2 = K\x2;
    end
    A = zeros(size(x1, 2), 9);
    for i=1:size(x1, 2)
        A(i, :) = ( kron(x1(:,i), x2(:,i)) )';
    end
    [U, S, V] = svd(A);
    
end
    %{
    buffer1 = ones(size(correspondences, 2));
    buffer0 = zeros(size(correspondences, 2));
    mask1 = [buffer1;buffer0;buffer1;buffer0];
    x1 = correspondences .* mask1;
    %}
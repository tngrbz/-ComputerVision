function sd = sampson_dist(F, x1_pixel, x2_pixel)
    % This function calculates the Sampson distance based on the fundamental matrix F
    sd = zeros(1, size(x1_pixel, 2));
    e3_hat = hat([0 0 1]');
    nom = sum(x2_pixel .* (F * x1_pixel), 1);
    denom1 = sum((e3_hat * F * x1_pixel).^2, 1);
    denom2 = sum((x2_pixel' * F * e3_hat).^2, 2);
    sd = nom.^2 ./ (denom1 + denom2');
end
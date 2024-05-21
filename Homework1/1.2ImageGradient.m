function [Fx, Fy] = sobel_xy(input_image)
    % In this function you have to implement a Sobel filter 
    % that calculates the image gradient in x- and y- direction of a grayscale image.
    disp(size(input_image))
    sobel_x = [1 0 -1;
        2 0 -2;
        1 0 -1];
    sobel_y = [1 2 1;
        0 0 0;
        -1 -2 -1];
    Fx = conv2(input_image, sobel_x, 'same');
    Fy = conv2(input_image, sobel_y, 'same');
    disp(size(Fx))
    disp(size(Fy))
end
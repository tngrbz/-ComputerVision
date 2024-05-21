function gray_image = rgb_to_gray(input_image)
    % This function is supposed to convert a RGB-image to a grayscale image.
    % If the image is already a grayscale image directly return it.
    if size(input_image,3) == 1
        gray_image = input_image;
        return
    end
    input_image = double(input_image);
    scaling = zeros(1,1,3);
    gray_image = zeros(size(input_image,1), size(input_image,2));
    scaling(1,1,:) = [0.299, 0.587, 0.114];
%     disp(size(input_image(1,1,:)));
%     disp(size(scaling));
    for i = 1:size(input_image,1)
        for j = 1:size(input_image,2)
            gray_image(i,j) = sum(input_image(i,j,:).*scaling);
        end
    end
    gray_image = uint8(gray_image);
end
%% Load images
Image1 = imread('sceneL.png');
IGray1 = rgb_to_gray(Image1);
imshow(IGray1);
Image2 = imread('sceneR.png');
IGray2 = rgb_to_gray(Image2);
imshow(IGray2);

%% Calculate Harris features
features1 = harris_detector(IGray1,'segment_length',9,'k',0.05,'min_dist',40,'N',50,'do_plot',false);
features2 = harris_detector(IGray2,'segment_length',9,'k',0.05,'min_dist',40,'N',50,'do_plot',false);

%% Correspondence estimation
correspondences = point_correspondence(IGray1,IGray2,features1,features2,'window_length',25,'min_corr',0.9,'do_plot',false);

%% Determine robust corresponding image points with the RANSAC algorithm
correspondences_robust = F_ransac(correspondences, 'tolerance', 0.04);

%% Visualize robust corresponding image points
figure
plot(correspondences_robust);

%% Calculate essential matrix
load('K.mat');
E = epa(correspondences, K);
disp(E);
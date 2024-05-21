classdef PointsReconstructor
    properties
        data_path                   % Path to the data
        camera_parameters           % Camera parameters (e.g., focal length, principal point)
        images                      % Cell array of images
        gui                         % gui to display progress
        view_set                    % View set
        view_set_iter               % Second view set that is refined iteratively
        points_3D                   % Reconstructed 3D points
        bad_connections             % Connections that should be deleted
        score_matrix                % Score matrix for the connections
        ref_image_id                % Index of the reference image
        second_ref_image            % Reference image id which should be used for the second triangulation (for the frames which have a big rotation angle)
        txt_camera_distances        % Camera distances between all cameras (given in the txt file)
        longest_path                % The longest path containing one connection from the reference
    end

    methods

        % Constructor
        function obj = PointsReconstructor(data_path, calibration_file_path, gui)
            obj.data_path = data_path;
            obj = obj.readCalibration(calibration_file_path);
            obj.gui = gui;
            obj.images = {};
            obj.view_set = imageviewset;
            obj.view_set_iter = imageviewset;
            obj.points_3D = [];
            obj.bad_connections = {};
            obj.score_matrix = sparse(numel(obj.images), numel(obj.images));
            obj.txt_camera_distances = [];
        end

        % Function to read the calibration from a file
        function obj = readCalibration(obj, file_path)
            % Open the file for reading
            fileID = fopen(file_path, 'r');
            % Read the file using textscan
            formatSpec = '%d %s %f %f %f %f %f %f';
            try
                data = textscan(fileID, formatSpec, 'CommentStyle', '#');
            catch ME
                fclose(fileID);
                error('Error reading the file: %s', ME.message);
            end

            % Close the file
            fclose(fileID);

            % Extract the camera information
            cameraID = data{1};
            model = data{2};
            width = data{3};
            height = data{4};
            params = [data{5}, data{6}, data{7}, data{8}];

            % Get the number of cameras
            numCameras = numel(cameraID);

            % Display the number of cameras
            fprintf('Number of cameras: %d\n\n', numCameras);

            % Display the camera information and create camera objects
            cameras = {};
            disp('Camera Information:');
            for i = 1:numCameras
                fprintf('Camera ID: %d\n', cameraID(i));
                fprintf('Model: %s\n', model{i});
                fprintf('Width: %.2f\n', width(i));
                fprintf('Height: %.2f\n', height(i));
                fprintf('Params: %.2f %.2f %.2f %.2f\n', params(i, :));
                fprintf('\n');
                % Create cameraParams object
                fx = params(i, 1);
                fy = params(i, 2);
                cx = params(i, 3);
                cy = params(i, 4);
                temp = cameraParameters('K', [fx 0 cx; 0 fy cy; 0 0 1], ...
                    'ImageSize', [height(i), width(i)]);
                cameras{end + 1} = temp;
            end

            % Add the parameters to the pointsReconstructer
            if length(cameras) == 1
                obj.camera_parameters = cameras{1};
            else
                disp(['More than one camera, please try another app!']);
            end
        end


        % Function to load the images
        function obj = loadImages(obj)
            
            % Create a file path pattern
            file_path_pattern = fullfile(obj.data_path, '*.jpg');
            % Get the list of jpg files
            file_list = dir(file_path_pattern);
            num_images = numel(file_list);

            % Loop over the images
            for i = 1:num_images
                % Get the file path
                file_name = file_list(i).name;
                file_path = fullfile(obj.data_path, file_name);
                % Load the image
                image = imread(file_path);
                % Add the image to the cell array
                obj.images{i} = im2gray(image);
            end
            fprintf('Found %d images\n', num_images);
        end

        % Function to extract the features and match them
        function obj = extractAndMatchFeatures(obj, k)
            
            % Loop over the images to extract the features
            for i = 1:numel(obj.images)
                % Get the image
                image = obj.images{i};
                % Detect SURF points and extract features for the current image
                points = detectSURFFeatures(image);
                [features, valid_points] = extractFeatures(image, points);
                % Add the view to the viewset object
                obj.view_set = addView(obj.view_set, i, 'Features', features, 'Points', valid_points);
            end

            % Loop over the images to match the features
            for i = 1:numel(obj.images)

                % Display progress
                fprintf('Matching features for image %d...\n', i);
                obj.gui.infoLabel.Text = ['Calculating coordinates...',newline,'Matching features for image ', num2str(i), '/', num2str(numel(obj.images))];
                drawnow
                % Initialize the number of matches for the current image with all the other images
                match_number = zeros(numel(obj.images), 1);
                % Initialize cell array to store the matched indices
                matched_indices = cell(numel(obj.images), 1);

                % Loop over the other images
                for j = 1:numel(obj.images)
                    if i == j
                        continue;
                    end

                    % Match the features
                    features1 = obj.view_set.Views{i, 3}{:}; % features are in the 3th column
                    features2 = obj.view_set.Views{j, 3}{:}; % features are in the 3th column
                    index_pairs = matchFeatures(features1, features2, 'Unique', true);
                    % Store the matched indices
                    matched_indices{j} = index_pairs;
                    % Update the number of matches
                    match_number(j) = size(index_pairs, 1);
                end

                % Get the indices of the k images with the most matches
                [~, indices] = maxk(match_number, k);

                % Create connections for the best matches
                for j = 1:numel(indices)
                    % Get the index of the other image
                    other_index = indices(j);
                    % Get the matched indices
                    index_pairs = matched_indices{other_index};
                    % Add the connection
                    if ~hasConnection(obj.view_set, i, other_index)
                        obj.view_set = addConnection(obj.view_set, i, other_index, 'Matches', index_pairs);
                    end

                end
            end
        end
        % Find the longest path in the dataset
        function obj = findLongestPath(obj)
            graph = createPoseGraph(obj.view_set);
            nodes = graph.Nodes.ViewId;
            optpath_length = 0;
            optpath = [];
            for i=1:1 % We always chose first kamera frame as reference
                v = dfsearch(graph,nodes(i),'edgetonew');
                if numel(v) == 0
                    continue
                end
                % Find unique elements
                [uniqueVec, ~, idx] = unique(v(:,1));
                counts = accumarray(idx, 1);
                % Find the repeating elements
                repeatingElements = uniqueVec(counts > 1);
                if numel(repeatingElements) == 0
                    cut_idx = length(v) + 1;
                else
                    cut_idx = length(v) + 1;
                    for c=1:numel(repeatingElements)
                        idx1 = find(v(:,1) == repeatingElements(c));
                        idx1 = idx1(2);
                        if cut_idx > idx1
                            cut_idx = idx1;
                        end
                    end
                end
                if cut_idx == 1
                    continue
                end
                path = [v(1:cut_idx-1,1); v(cut_idx-1, 2)]';
                path_length = length(path);
                if path_length > optpath_length
                    optpath = path;
                    optpath_length = path_length;
                end
            end
            obj.longest_path = optpath;
            
        end


        % Function to link the cameras together
        function obj = findRelativeCameraPoses(obj)
            number_of_connections = size(obj.view_set.Connections,1);
            
            for i=1:number_of_connections
                viewId1 = obj.view_set.Connections{i,1} ;
                viewId2 = obj.view_set.Connections{i,2} ; 

                validPoints1 = obj.view_set.Views{viewId1,4}{1} ;
                validPoints2 = obj.view_set.Views{viewId2,4}{1} ; 
                
                indexPairs = obj.view_set.Connections{i,5}{1} ;
                % MatchedPoints contain the indeces of the matched features in their
                % corresponding image.
                matchedPoints1 = validPoints1(indexPairs(:, 1), :);
                matchedPoints2 = validPoints2(indexPairs(:, 2), :);
                
                ransacIterations = 10000;  % Number of iterations for RANSAC
                ransacInlierDistanceThreshold = 0.1;  % Maximum allowed distance between a point and its corresponding epipolar line for inliers

                % Perform RANSAC to estimate the fundamental matrix
                [E, inlierIdx, status] = estimateEssentialMatrix(matchedPoints1, matchedPoints2, obj.camera_parameters, ...
                    "MaxNumTrials", ransacIterations, 'MaxDistance', ransacInlierDistanceThreshold);

                % Filter the matched points using the RANSAC inlier indices
                % You further eliminate the features which are more distant from each other
                % than threshold.
                inlierPoints1 = matchedPoints1(inlierIdx, :);
                inlierPoints2 = matchedPoints2(inlierIdx, :);

                filtered_matches = [indexPairs(inlierIdx, 1) indexPairs(inlierIdx, 2)];
                obj.view_set = updateConnection(obj.view_set, viewId1, viewId2, 'Matches', filtered_matches);

                % Recover the relative camera pose (rotation and translation) from the essential matrix
                intrinsics = obj.camera_parameters.Intrinsics ;
                try
                    [relativePose,validPointsFraction] = estrelpose(E,intrinsics,inlierPoints1,inlierPoints2);
                    % Check if validPointsFraction is not too small, otherwise delete the connection
                    % (see documentation of estrelpose for the meaning of validPointsFraction)
                    if validPointsFraction < 0.9 || status == 1 || status == 2
                        obj.bad_connections{end + 1} = [viewId1,viewId2];
                        continue;
                    end

                    % Add the validPointsFraction to the score matrix
                    obj.score_matrix(viewId1,viewId2) = validPointsFraction;
    
                    % Multiply the translation vector by the distance between the cameras (given in the txt file)
                    %relativePose.Translation = relativePose.Translation * obj.txt_camera_distances(viewId1, viewId2);
                
                    obj.view_set = updateConnection(obj.view_set,viewId1,viewId2,relativePose);
                catch
                    obj.bad_connections{end + 1} = [viewId1,viewId2];
                end
                
            end

            % Delete the connections that have to be deleted
            for i=1:numel(obj.bad_connections)
                obj.view_set = deleteConnection(obj.view_set,obj.bad_connections{i}(1),obj.bad_connections{i}(2));
            end

            % Delete the views that have no connections or just one connection
            temp = obj.view_set.Connections.ViewId1;
            for i=1:obj.view_set.NumViews
                if size(temp(temp == i),1) <= 1
                    %obj.view_set = deleteView(obj.view_set,i);
                end
            end
        
        end

        % Calculate absolute poses with respect to the first image.
        function obj = get_absolute_poses(obj)
            Im1 = obj.view_set.Views.ViewId(1,1);
            % Created a queue object to iterate over comparisons
            % It contains the indices of Connections, where initially the
            % first image is compared to. Afterwards you add other
            % comparisons for other images.
            queue = find(obj.view_set.Connections.ViewId1(:) == Im1);
            while true
                if isempty(queue)
                    break
                end
                idx_comp = queue(1); % Comparison index
                queue = queue(2:end); % then I pop the first element of the queue.
                Im1 = obj.view_set.Connections.ViewId1(idx_comp,1); % Im1 is the first index of comparison.
                Im2 = obj.view_set.Connections.ViewId2(idx_comp,1); % Im2 is the second index of the comparison.
                % Get tabular indices (as logical vectors)
                tab_idx1 = obj.view_set.Views.ViewId == Im1;
                tab_idx2 = obj.view_set.Views.ViewId == Im2;
                A_0 = [eye(3), zeros(3,1); zeros(1,3), 1]; % This is the initial transmission

                if ~all(obj.view_set.Views.AbsolutePose(tab_idx1).A(:) == A_0(:)) && ~all(obj.view_set.Views.AbsolutePose(tab_idx2).A(:) == A_0(:))
                    % If both images have an absolute position, we skip
                    continue
                elseif Im2 == obj.view_set.Views.ViewId(1,1) && ~all(obj.view_set.Views.AbsolutePose(tab_idx1).A(:) == A_0(:))
                    % If Im2 is the reference image, and Im1 has an
                    % absolute value, we skip
                    continue
                end

                A = obj.view_set.Connections.RelativePose{idx_comp}.A; % Relative Transformation
                % One thing to notice is that Relative Pose is a cell and
                % Absolute Pose is a 4x1 datatype
                abs_pose = obj.view_set.Views.AbsolutePose(tab_idx1).A; % Absolute Position of the first image.
                abs_pose2 = rigidtform3d(A*abs_pose); % Multiply from the left side to get abs_pose of the secondary image
                obj.view_set = updateView(obj.view_set,Im2,abs_pose2); % Update the view_set again.

                buffer = find(obj.view_set.Connections.ViewId1(:) == Im2);
                k = 1;
                for i=1:length(buffer)
                    if any(queue == buffer(k))
                        buffer(k) = [];
                        continue
                    end
                    k = k+1;
                end
                queue = [queue; buffer];
            end
        end

        % Function to get the absolute poses automatically
        function obj = optimizeAbsolutePoses(obj)
            obj.view_set = optimizePoses(obj.view_set, 'MaxIterations', 5000, 'Tolerance', 1e-10);
        end

        % Function to plot the camera poses of given views
        function obj = plot_camera_poses(obj, view_set)

            % Get camera poses
            camera_poses = poses(view_set);

            % Plot the camera poses
            figure;
            plotCamera(camera_poses, 'Size', 0.1);
            grid on
            campos([5, -10, -30]);
            camup([0, -1, 0]);
            % Label axes
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
    
        end

        %% Function to plot camera poses and 3D points in one figure
        function obj = plot_camera_poses_and_points(obj, view_set, viewIds, threshold)

            % Get camera poses
            camera_poses = poses(view_set, viewIds);

            % Plot the camera poses
            figure;
            plotCamera(camera_poses, 'Size', 0.1); hold on;
            grid on
            campos([5, -10, -30]);
            camup([0, -1, 0]);

            % Remove all points that do not have coordinates between -30 and 30
            % (outliers)
            obj.points_3D(obj.points_3D(:,1) < -threshold | obj.points_3D(:,1) > threshold, :) = [];
            obj.points_3D(obj.points_3D(:,2) < -threshold | obj.points_3D(:,2) > threshold, :) = [];
            obj.points_3D(obj.points_3D(:,3) < -threshold | obj.points_3D(:,3) > threshold, :) = [];

            % % Switch the y and z coordinates of the points
            % temp = obj.points_3D(:,2);
            % obj.points_3D(:,2) = obj.points_3D(:,3);
            % obj.points_3D(:,3) = temp;
            obj.points_3D = pointCloud(obj.points_3D);
            obj.points_3D = pcdenoise(obj.points_3D, "Threshold", 0.1);
            pcshow(obj.points_3D);
            camup([0, -1, 0]);
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title('Camera poses and 3D Points');
            hold off;
        end

        % Function to get an initial estimate of the points
        function obj = get_initial_estimate(obj)
            point_tracks = findTracks(obj.view_set);
            camera_poses = poses(obj.view_set);
            intrinsics = obj.camera_parameters.Intrinsics;
            obj.points_3D = triangulateMultiview(point_tracks, camera_poses, intrinsics);
            % disp(obj.points_3D)
            % disp(size(obj.points_3D))
            % disp(class(obj.points_3D))
        end

        % Function to refine the points
        function obj = refine_points(obj)
            point_tracks = findTracks(obj.view_set);
            [points_3D_refinded, refined_poses] = bundleAdjustment(obj.points_3D, point_tracks, poses(obj.view_set), obj.camera_parameters.Intrinsics);
            obj.points_3D = points_3D_refinded;
       

            % Plot the camera poses
            figure;
            plotCamera(refined_poses, 'Size', 0.1, 'Color', [0 0 1], 'Opacity', 0.05); hold on;
            grid on
            campos([5, -10, -30]);
            camup([0, -1, 0]);
        
            %pcshow(obj.points_3D);
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title('Camera poses and 3D Points');
            hold off;
        end

        %% Function to get 3D points from a given imageviewset
        function obj = get_3D_points_from_imageviewset(obj, viewIds)

            [xyzPoints, vSet] = obj.get_3D_points_iteratively(viewIds);
            
            % Add the 3D points and the viewset to the object
            obj.points_3D = xyzPoints;
            obj.view_set_iter = vSet;
        end

        % Function that creates a viewset iteratively from a given path and returns 3D points
        function [xyzPoints, vSet] = get_3D_points_iteratively(obj, path)
            
            % Create a viewset
            vSet = imageviewset;

            % Add the reference image to the viewset
            ref_id = path(1);
            ref_view = findView(obj.view_set, ref_id);
            vSet = addView(vSet, path(1), rigidtform3d, ...
                            'Points', ref_view.Points{:}, 'Features', ref_view.Features{:}');
            
            % Loop over the path
            for i = 1:numel(path) - 1

                % Check if the view is already in the viewset
                if hasView(vSet, path(i+1))
                    continue;
                end
                
                % Get absolute pose of the previous image
                previous_pose = poses(vSet, path(i)).AbsolutePose;
                % Get relative pose from the previous image to the current image
                relative_pose = findConnection(obj.view_set, path(i), path(i+1)).RelativePose{:};
                % Compute the absolute pose of the current image
                current_pose = rigidtform3d(previous_pose.A * relative_pose.A);

                % Update the viewset
                current_view = findView(obj.view_set, path(i+1));
                vSet = addView(vSet, path(i+1), current_pose, ...
                                'Points', current_view.Points{:}, 'Features', current_view.Features{:}');
                current_connection = findConnection(obj.view_set, path(i), path(i+1));
                vSet = addConnection(vSet, path(i), path(i+1), relative_pose, 'Matches', current_connection.Matches{:});
                
                % Find arguments for triangulateMultiview
                points_tracks = findTracks(vSet);
                camera_poses = poses(vSet);
                intrinsics = obj.camera_parameters.Intrinsics;

                % Triangulate the points
                [xyzPoints, errors, valid_idx] = triangulateMultiview(points_tracks, camera_poses, intrinsics);

                % Do bundle adjustment
                [xyzPoints, camera_poses, reprojection_errors] = bundleAdjustment(xyzPoints, points_tracks, ...
                                                                    camera_poses, intrinsics, PointsUndistorted=true);
                % Remove all points that are behind the cameras
                % Remove all points that have big reprojection errors
                reprojection_errors = reprojection_errors < 10;
                filter = reprojection_errors & valid_idx;
                xyzPoints = xyzPoints(filter, :);
                % Update the camera poses
                vSet = updateView(vSet, camera_poses);                                    
            end
        end

        % Function to plot the digraph of the viewset
        function plot_graph(obj)
            figure;
            plot(obj.view_set, ShowViewIds='on');
        end

        % Function to plot the shortest paths from a reference image to all other images
        function plot_shortest_paths(obj, ref_image_id)
            figure;
            graph = createPoseGraph(obj.view_set);
            tree = shortestpathtree(graph, ref_image_id, "Method", "unweighted");
            p = plot(graph);
            highlight(p, tree, "EdgeColor", "r", "LineWidth", 3);
        end

        % Function to get the 3D points iteratively from a given imageviewset
        function xyzPoints = get_3D_points(obj)

            % Load the images
            disp('Loading images...');
            obj.gui.infoLabel.Text = ['Calculating coordinates...',newline,'Loading images'];
            drawnow
            obj = obj.loadImages();

            % Extract the features
            disp('Extracting features...');
            obj.gui.infoLabel.Text = ['Calculating coordinates...',newline,'Extracting features'];
            drawnow
            obj = obj.extractAndMatchFeatures(3);

            % Estimate relative poses
            disp('Estimating relative poses...');
            obj.gui.infoLabel.Text = ['Calculating coordinates...',newline,'Estimating relative poses'];
            drawnow
            obj = obj.findRelativeCameraPoses();

            % Get the points iteratively from the longest path (hardcoded for now)
            disp('Getting 3D points...');
            obj.gui.infoLabel.Text = ['Calculating coordinates...',newline,'Getting 3D-points'];
            drawnow
            obj = findLongestPath(obj);
            path = obj.longest_path;
            obj = obj.get_3D_points_from_imageviewset(path);
            %obj = obj.svdCorrection();
            % Return the 3D points
            xyzPoints = obj.points_3D;
        end 

       function obj = CalculateCameraDistances(obj, file_path)
            % Read the contents of the file
            data = readlines(file_path);

            %initialize array which contains the ImageID and the rigidt3Dtransform
            rigidt_struct = struct();

            for i = 1:size(data, 1)
    
                %split each line of strings to individual strings
                splitelements = split(data(i));
    
                %convert all the element to type double
                array = str2double(splitelements(1,1)); 
    
                %filter out all lines that dont contain valuable information out
                naturalNumber = array(array == floor(array));
    
                %access this line
                if naturalNumber > 0
                    %extract the elements 
                    TX = str2double(splitelements(6,1));
                    TY = str2double(splitelements(7,1));
                    TZ = str2double(splitelements(8,1));
        
                    %calculate translation vector
                    vec_t = [TX TY TZ];
         
                    %Store all the ImageId's and Rigidt Objects in the struct rigidt_struct
                    fieldName = sprintf('RigidtObject%d', naturalNumber);
                    rigidt_struct.(fieldName).ImageID = naturalNumber;
                    rigidt_struct.(fieldName).Translation = vec_t; 

                end 
            end

            %Iterate over the other images 
            for i = 1:length(fieldnames(rigidt_struct)) 
            
                fieldName_i = sprintf('RigidtObject%d', i);
                
                %Get T from struct 
                T_i = rigidt_struct.(fieldName_i).Translation; 
                
                for j = 1:length(fieldnames(rigidt_struct)) 
                    if j == i
                        break
                    else
                        fieldName_j = sprintf('RigidtObject%d', j);

                        %Get T from struct  
                        T_j = rigidt_struct.(fieldName_j).Translation; 
               
                        %Calculate vector starting from T_i and pointing to T_j
                        T_res = T_j - T_i;

                        %Calculate the Norm 
                        norm_T_res = norm(T_res);
            
                        %Store into struct 
                        obj.txt_camera_distances(i, j) =norm_T_res;
                        obj.txt_camera_distances(j, i) =norm_T_res;
                    end 
                end
            end
       end  
       function obj = svdCorrection(obj)
               X = obj.points_3D';
               [U S V] = svd(X);
               if (det(U) ~= 1)
                    U = U*[1 0 0;0 1 0;0 0 -1];
               end
               obj.points_3D = U'*X;
               obj.points_3D = obj.points_3D';
               temp = obj.points_3D(:,2);
               obj.points_3D(:,2) =obj.points_3D(:,3);
               obj.points_3D(:,3) = temp ;
       end


       

    end
end 
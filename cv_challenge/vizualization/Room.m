classdef Room
    properties
        Boxes
        gui
        fullPointCloud
    end
    
    methods
        function obj = Room(points, app)
            
            obj.gui = app;

            % Calculate the Z-scores
            zScores = zscore(points);
            
            % Identify outliers
            outliers = any(abs(zScores) > 2.5, 2);
            
            % Remove outliers
            pointsArray = points(~outliers, :);
            size(pointsArray)
            
            obj.showWholePointcloud(pointCloud(pointsArray))
            tmpPointCloud = pcdenoise(pointCloud(pointsArray));
            obj.showWholePointcloud(tmpPointCloud)
            vectorArray = num2cell(tmpPointCloud.Location, 2);
            
            % Convert the cell array of 1x3 vectors back to a matrix
            pointsArray = cell2mat(vectorArray);
            size(pointsArray)
            maxDistance = 0.1;
            
            % Plane fitting and segmentation using RANSAC
            [model1,inlierIndices,outlierIndices] = pcfitplane(pointCloud(pointsArray),maxDistance);
            inlierPoints = pointsArray(inlierIndices,:);
            outlierPoints = pointsArray(outlierIndices,:);
            
            % Initial number of clusters
            k = 30;  
            
            % Do k-means on the inlier points
            [idx, ~] = kmeans(inlierPoints, k);
            
            % Calculate the silhouette score for the inlier points
            s = silhouette(inlierPoints, idx);
            avgSilhouette = mean(s);
            
            bestK = k;

            % Find the optimal number of clusters for the inlier points
            while true && k < 120
                k = k + 1;
                [idx, ~] = kmeans(inlierPoints, k);
                s = silhouette(inlierPoints, idx);
                newAvgSilhouette = mean(s);
                
                if newAvgSilhouette < avgSilhouette
                    % The silhouette score has decreased, so the previous number of clusters was better
                    
                else
                    avgSilhouette = newAvgSilhouette;
                    bestK = k;
                end
            end
            display(bestK)
            display(avgSilhouette)
            % Do k-means on the inlier points
            [idx, ~] = kmeans(inlierPoints, bestK);
            
            outlierK = max(1, round(sqrt(size(outlierPoints, 1))));
            if size(outlierPoints, 1) > 0
                % Do the same for the outlier points but with a fixed number of clusters (you could also optimize this separately)
                [idxOutlier, ~] = kmeans(outlierPoints, outlierK);
            else
                outlierK = 0;
                idxOutlier = [];
            end

            % Create a cell array to store the pointClouds
            inlierPointClouds = cell(1, bestK);
            outlierPointClouds = cell(1, outlierK);
            
            % Create a pointCloud for each inlier cluster and filter out the ceiling
            for i = 1:bestK
                tempPointCloud = pointCloud(inlierPoints(idx == i, :));
                % Check if the mean Z value of the point cloud is less than -1
                % if mean(tempPointCloud.Location(:, 3)) < -2
                    % Calculate Z-scores for the point cloud
                    tempZScores = zscore(tempPointCloud.Location);
                    % Identify outliers in the point cloud
                    tempOutliers = any(abs(tempZScores) > 2.5, 2);
                    % Remove outliers from the point cloud
                    tempPointCloud = pointCloud(tempPointCloud.Location(~tempOutliers, :));
                    inlierPointClouds{i} = tempPointCloud;
                % end
                %inlierPointClouds{i} = tempPointCloud;
            end
            
            % Create a pointCloud for each outlier cluster and filter out the ceiling
            for i = 1:outlierK
                tempPointCloud = pointCloud(outlierPoints(idxOutlier == i, :));
                
                % Check if the mean Z value of the point cloud is less than -1
                % if mean(tempPointCloud.Location(:, 3)) < -2
                    % Calculate Z-scores for the point cloud
                    tempZScores = zscore(tempPointCloud.Location);
                    % Identify outliers in the point cloud
                    tempOutliers = any(abs(tempZScores) > 2.5, 2);
                    % Remove outliers from the point cloud
                    tempPointCloud = pointCloud(tempPointCloud.Location(~tempOutliers, :));
                    outlierPointClouds{i} = tempPointCloud;
                % end
                outlierPointClouds{i} = tempPointCloud;
            end
            
            % Combine inlier and outlier point clouds

            % Filter out empty cells from inlierPointClouds and outlierPointClouds
            inlierPointClouds = inlierPointClouds(~cellfun('isempty', inlierPointClouds));
            outlierPointClouds = outlierPointClouds(~cellfun('isempty', outlierPointClouds));
            
            pointClouds = [inlierPointClouds, outlierPointClouds];
            
            % Create Box objects for each cluster and add them to the Room
            obj.Boxes = cell(1, numel(pointClouds));
            for i = 1:numel(pointClouds)
                try
                    obj.Boxes{i} = Box(pointClouds{i});
                catch exception
                    fprintf('Error while creating box for point cloud %d:\n', i);
                    disp(exception);
                end
            end
            
            % Merge point clouds using a for loop
            mergedPoints = [];
            mergedColor = [];
            mergedIntensity = [];
            
            for i = 1:numel(pointClouds)
                mergedPoints = [mergedPoints; pointClouds{i}.Location];
                mergedColor = [mergedColor; pointClouds{i}.Color];
                mergedIntensity = [mergedIntensity; pointClouds{i}.Intensity];
            end
            
            % Create the merged point cloud object
            mergedPC = pointCloud(mergedPoints, 'Color', mergedColor, 'Intensity', mergedIntensity);
            obj.fullPointCloud = mergedPC;

        end

        function showWholePointcloud(obj, pointcloud)
            axes = obj.gui.RoomFigure;
            hold(axes, 'on');
            color = rand(1,3); % Generate a random RGB color
            pcshow(pointcloud.Location, color, 'Parent', axes);
        end

        function showWalls(obj)
            % Retrieve all points
            points = obj.fullPointCloud.Location;
            axes = obj.gui.RoomFigure;
            
            boundaryIndices = boundary(points(:,1), points(:,2), 0.3);
            outerWallPoints = points(boundaryIndices,:);
            
            % The wall's height (Z direction) is from min Z to max Z of all points
            minZ = min(outerWallPoints(:, 3));
            maxZ = max(outerWallPoints(:, 3));
            
            hold(axes, 'on');
            % Create the walls
            for i = 1:size(outerWallPoints, 1)
                % create two polygons for each pair of points: one at the floor (minZ) and one at the ceiling (maxZ)
                
                j = mod(i, size(outerWallPoints, 1)) + 1;  % This will make j loop back to 1 when i reaches the end

                %if pdist2(outerWallPoints(i, 1:2), outerWallPoints(j, 1:2)) > 9
                %    continue; % Skip the current iteration if the distance is larger than 5
                %end
            
                fill3([outerWallPoints(i, 1), outerWallPoints(j, 1), outerWallPoints(j, 1), outerWallPoints(i, 1)], ...
                      [outerWallPoints(i, 2), outerWallPoints(j, 2), outerWallPoints(j, 2), outerWallPoints(i, 2)], ...
                      [minZ, minZ, maxZ, maxZ], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'k', 'Parent', axes);
            end
        end

        function showWallPlot(obj)
            tmpPoints = obj.fullPointCloud.Location;
            points = [tmpPoints(:,1), tmpPoints(:,2)];

            boundary_indexs = boundary(points(:,1), points(:,2), 0.4);
            boundary_points = points(boundary_indexs,:);


            boundary_indexs2 = boundary(points(:,1), points(:,2), 0.2);
            boundary_points2 = points(boundary_indexs2,:);
            
            % Plot the original points and the simplified polygon
            figure;
            plot(points(:, 1), points(:, 2), 'bo');
            hold on;
            %plot(simplified_points(:, 1), simplified_points(:, 2), 'ro-');
            plot(boundary_points(:, 1), boundary_points(:, 2), 'ro-');
            plot(boundary_points2(:, 1), boundary_points2(:, 2), 'ko-');
            
            legend('Original Points','boundary 0,4', 'boundary 0.2');
        end
        
        function show(obj)
            axes = obj.gui.RoomFigure;
            hold(axes, 'on');
            for i = 1:numel(obj.Boxes)
                %obj.Boxes{i}.showPoints(axes);
                obj.Boxes{i}.showBox(axes, @(src, event) obj.boxClicked(src, event, i));  % Pass the box index to the callback function
            end
            %{
            % Only for testing:
            obj.Boxes{1}.showPoints(axes);
            obj.Boxes{1}.showBox(axes, @(src, event) obj.boxClicked(src, event, 1, axes));
            
            obj.Boxes{2}.showPoints(axes);
            obj.Boxes{2}.showBox(axes, @(src, event) obj.boxClicked(src, event, 2, axes));
            
            obj.Boxes{3}.showPoints(axes);
            obj.Boxes{3}.showBox(axes, @(src, event) obj.boxClicked(src, event, 3, axes));
            obj.Boxes{4}.showPoints(axes);
            obj.Boxes{4}.showBox(axes, @(src, event) obj.boxClicked(src, event, 4, axes));
            obj.Boxes{5}.showBox(axes, @(src, event) obj.boxClicked(src, event, 5, axes));
            obj.Boxes{6}.showBox(axes, @(src, event) obj.boxClicked(src, event, 6, axes));
            obj.Boxes{7}.showBox(axes, @(src, event) obj.boxClicked(src, event, 7, axes));
            %}
            view(axes, 3);
            axis(axes, 'equal');
            grid(axes, 'on');
            hold(axes, 'off');
        end
        
        function boxClicked(obj, src, event, index)
            axes = obj.gui.RoomFigure;
            label1 = obj.gui.Box1SizeLabel;
            label2 = obj.gui.Box2SizeLabel;

            persistent selectedBox selectedFace lineHandle textHandle

            if ~isempty(lineHandle) && isvalid(lineHandle)
                delete(lineHandle);
                lineHandle = [];
            end

            if ~isempty(textHandle) && isvalid(textHandle)
                delete(textHandle);
                textHandle = [];
            end
            
            display(selectedBox)

            if isempty(selectedBox)
                selectedBox = obj.Boxes{index};
                selectedFace = str2double(src.Tag);  % Find the clicked face based on the tag
                 % Display the size of the first box on the label
                label1.Text = selectedBox.getStringSize();
                label2.Text = 'Select second Box';
            else
                box1 = selectedBox;
                face1 = selectedFace;
                box2 = obj.Boxes{index};
                face2 = str2double(src.Tag);
                
                display(face1)
                display(face2)

                box1.FaceCenters
                box2.FaceCenters

                point1 = box1.FaceCenters(face1, :);
                point2 = box2.FaceCenters(face2, :);

                distance = norm(point1 - point2);

                fprintf('The distance between the selected boxes is %.2f\n', distance);
                % Display the size of the first box on the label
                label2.Text = [box2.getStringSize(), "m"];

                
                midpoint = (point1 + point2) / 2;
                lineHandle = line(axes, [point1(1), point2(1)], [point1(2), point2(2)], [point1(3), point2(3)], 'Color', 'r');
                textHandle = text(axes, midpoint(1), midpoint(2), midpoint(3), sprintf('%.2f', distance), 'Color', 'r');

                selectedBox = [];
            end
        end
    end
end
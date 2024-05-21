classdef Box
    properties
        MinCorner
        MaxCorner
        Center
        PointCloud
        FaceCenters
    end

    properties (Dependent = true)
        Dimensions
        Volume
        StringSize
    end
    
    methods
        function obj = Box(pointCloud)
            obj.PointCloud = pointCloud;
            obj.MinCorner = min(pointCloud.Location, [], 1);
            obj.MaxCorner = max(pointCloud.Location, [], 1);
            obj.Center = mean([obj.MinCorner; obj.MaxCorner]);

            % Calculate the box vertices
            V = [obj.MinCorner(1), obj.MinCorner(2), obj.MinCorner(3);
                 obj.MaxCorner(1), obj.MinCorner(2), obj.MinCorner(3);
                 obj.MaxCorner(1), obj.MaxCorner(2), obj.MinCorner(3);
                 obj.MinCorner(1), obj.MaxCorner(2), obj.MinCorner(3);
                 obj.MinCorner(1), obj.MinCorner(2), obj.MaxCorner(3);
                 obj.MaxCorner(1), obj.MinCorner(2), obj.MaxCorner(3);
                 obj.MaxCorner(1), obj.MaxCorner(2), obj.MaxCorner(3);
                 obj.MinCorner(1), obj.MaxCorner(2), obj.MaxCorner(3)];
            
            % Calculate face centers
            obj.FaceCenters = [(V(1, :) + V(2, :) + V(3, :) + V(4, :)) / 4;
                               (V(5, :) + V(6, :) + V(7, :) + V(8, :)) / 4;
                               (V(1, :) + V(2, :) + V(6, :) + V(5, :)) / 4;
                               (V(3, :) + V(4, :) + V(8, :) + V(7, :)) / 4;
                               (V(1, :) + V(4, :) + V(8, :) + V(5, :)) / 4;
                               (V(2, :) + V(3, :) + V(7, :) + V(6, :)) / 4];

        end

        function dimensions = getDimensions(obj)
            dimensions = abs(obj.MaxCorner - obj.MinCorner);
        end

        function stringSize = getStringSize(obj)
            dimensions = abs(obj.MaxCorner - obj.MinCorner);
            display(dimensions)
            stringSize = num2str(dimensions(1),2) + "m x" + num2str(dimensions(2),2) + "m x" + num2str(dimensions(3),2) + "m";
            %stringSize = sprintf('%2g x %2g x %2g', dimensions(1), dimensions(2), dimensions(3)) ;
        end
        
        function vol = getVolume(obj)
            dimensions = abs(obj.MaxCorner - obj.MinCorner);
            vol = dimensions(1) * dimensions(2) * dimensions(3);
        end

        
        function showBox(obj, axes, callbackFcn)
            V = [obj.MinCorner(1), obj.MinCorner(2), obj.MinCorner(3);
                 obj.MaxCorner(1), obj.MinCorner(2), obj.MinCorner(3);
                 obj.MaxCorner(1), obj.MaxCorner(2), obj.MinCorner(3);
                 obj.MinCorner(1), obj.MaxCorner(2), obj.MinCorner(3);
                 obj.MinCorner(1), obj.MinCorner(2), obj.MaxCorner(3);
                 obj.MaxCorner(1), obj.MinCorner(2), obj.MaxCorner(3);
                 obj.MaxCorner(1), obj.MaxCorner(2), obj.MaxCorner(3);
                 obj.MinCorner(1), obj.MaxCorner(2), obj.MaxCorner(3)];
        
            F = [1, 2, 3, 4; % Bottom face
                 5, 6, 7, 8; % Top face
                 1, 2, 6, 5; % Front face
                 3, 4, 8, 7; % Back face
                 1, 4, 8, 5; % Left face
                 2, 3, 7, 6]; % Right face
        
            for i = 1:size(F, 1)
                faceHandle = fill3(axes, V(F(i, :), 1), V(F(i, :), 2), V(F(i, :), 3), 'g', 'FaceAlpha', 0.1, 'ButtonDownFcn', callbackFcn); 
                faceHandle.Tag = num2str(i);  % Add a unique tag to each face
            end
        end

        function showPoints(obj, axes)
            color = rand(1,3); % Generate a random RGB color
            pcshow(obj.PointCloud.Location, color, 'Parent', axes);
        end

        function points = surfaceGrid(obj, numPoints)
            [x, y, z] = meshgrid(linspace(obj.MinCorner(1), obj.MaxCorner(1), numPoints), ...
                                 linspace(obj.MinCorner(2), obj.MaxCorner(2), numPoints), ...
                                 linspace(obj.MinCorner(3), obj.MaxCorner(3), numPoints));
            points = [x(:), y(:), z(:)];
            onXSurface = ismember(points(:, 1), [obj.MinCorner(1), obj.MaxCorner(1)]);
            onYSurface = ismember(points(:, 2), [obj.MinCorner(2), obj.MaxCorner(2)]);
            onZSurface = ismember(points(:, 3), [obj.MinCorner(3), obj.MaxCorner(3)]);
            points = points(onXSurface | onYSurface | onZSurface, :);
        end
    end
end
%}
function [samplingPositions, gridDimension] = generateMeasurementGrid(minLat, maxLat, minLon, maxLon, gridSpacing, ueHeight, simpleBuildings)
% GENERATEMEASUREMENTGRID Generates grid strictly within Min/Max bounds.

% Show whats happening
fprintf('Generating Grid from Lat: %.4f-%.4f, Lon: %.4f-%.4f\n', minLat, maxLat, minLon, maxLon);

% Earth Circumference approx
R = 6378137;
circLat = 2 * pi * R; % ~40,075km

% Calculate Step Sizes in Degrees
% Latitude: Constant (~111km per degree)
degPerMeterLat = 360 / circLat;

stepLat = gridSpacing * degPerMeterLat;

% Longitude: Depends on Latitude (shrinks North direction)
% Average latitude to calculate the longitude scaling
avgLat = (minLat + maxLat) / 2;
radiusAtLat = R * cosd(avgLat);
circLon = 2 * pi * radiusAtLat;
degPerMeterLon = 360 / circLon;

stepLon = gridSpacing * degPerMeterLon;

% Generate Vectors
latVec = minLat:stepLat:maxLat;
lonVec = minLon:stepLon:maxLon;

[lonGrid, latGrid] = meshgrid(lonVec, latVec); % Combine lat/lon vectors to grid

% Initialize Mask
validMask = true(size(latGrid));

% Filter Buildings ot of samplePoints
if nargin > 6 && ~isempty(simpleBuildings) % nargin = #function imputs, check if all are given
    fprintf('Processing exclusion zones... ');

    % Flatten grid for analysis
    qLat = latGrid(:);
    qLon = lonGrid(:);

    % Convert struct to polyshape vector
    numBuildings = length(simpleBuildings);
    polyVec = repmat(polyshape, numBuildings, 1);

    for i = 1:numBuildings
        % [Lon, Lat] corresponds to [X, Y]
        % Create rects for buildings with polyshape
        polyVec(i) = polyshape(simpleBuildings(i).Lon,...
            simpleBuildings(i).Lat,...
            'KeepCollinearPoints', true);
    end

    % Merge all buildings together
    allBuildingsPoly = union(polyVec);

    % Check Points
    inPoly = isinterior(allBuildingsPoly, qLon, qLat);
    validMask(inPoly) = false; % Set points which are in buildings/roofs not good

    fprintf('Removed %d points.\n', sum(inPoly));
end

% Final Formatting
% Just keep valid points (not in buildings)
latFlat = latGrid(validMask);
lonFlat = lonGrid(validMask);
hFlat = ones(size(latFlat)) * ueHeight; % Ue always 1.5m high, works with rxsite?

samplingPositions = [latFlat, lonFlat, hFlat];

% Return dimension info if needed for reshaping later
% and building the plot
gridDimension.rows = size(latGrid, 1); 
gridDimension.cols = size(latGrid, 2); 
gridDimension.latGrid = latGrid; % Needed for plotting axis
gridDimension.lonGrid = lonGrid; % Needed for plotting axis
gridDimension.validMask = validMask; % CRITICAL for reconstruction
gridDimension.validMask = validMask; % for reconstruction
gridDimension.minLat = minLat;
gridDimension.minLon = minLon;
gridDimension.stepLat = stepLat;
gridDimension.stepLon = stepLon;

fprintf('Grid ready: %d points.\n', length(samplingPositions));
end
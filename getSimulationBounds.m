function [minLat, maxLat, minLon, maxLon, gridSpacing, buildingStruct] = getSimulationBounds(measureRx, osmFile, centerLat, centerLon, lambda)
% GETSIMULATIONBOUNDS Determines the simulation area boundaries
%
% Inputs:
%   measureRx  - 0 for Full Map (OSM Bounds), 1 for 1x1m Local Box
%   osmFile    - String path to the .osm file (used for Full Map bounds)
%   centerLat  - Receiver Latitude (used for Local Box center)
%   centerLon  - Receiver Longitude (used for Local Box center)
%   lambda     - Wavelength (used for grid spacing in Local Box)
%
% Outputs:
%   minLat, maxLat, minLon, maxLon - The boundaries of the simulation
%   gridSpacing - The distance between points (m)

if measureRx == 0
    % Full Map Overview
    fprintf('Mode: Full Map Overview selected...\n');

    gridSpacing = 2; % Try 1 meter like Tam (averaged) maybe with Lambda/2 directly around UE (m)

    % Read Bounds directly from the XML Header
    xDoc = xmlread(osmFile);
    boundsItems = xDoc.getElementsByTagName('bounds');
    boundsTag = boundsItems.item(0);

    minLat = str2double(boundsTag.getAttribute('minlat'));
    minLon = str2double(boundsTag.getAttribute('minlon'));
    maxLat = str2double(boundsTag.getAttribute('maxlat'));
    maxLon = str2double(boundsTag.getAttribute('maxlon'));

    fprintf('Using OSM File Bounds:\n Lat: %.4f - %.4f\n Lon: %.4f - %.4f\n', ...
        minLat, maxLat, minLon, maxLon);


else
    % High Res 1x1m Square around Center RX
    % Should show fading effects
    fprintf('Mode: High Resolution Local Sampling (Lambda/2) selected.\n');

    gridSpacing = lambda / 10; % High resolution spacing

    % Geodesy parameters
    R = 6378137; % Earth Radius in meters
    boxSize = 1; % 1 meter box
    halfSide = boxSize / 2;

    % Calculate Latitude Offset
    deltaLat = (halfSide / R) * (180 / pi);

    % Calculate Longitude Offset
    deltaLon = (halfSide / (R * cosd(centerLat))) * (180 / pi);

    % Set Bounds centered on the provided coordinates
    minLat = centerLat - deltaLat;
    maxLat = centerLat + deltaLat;
    minLon = centerLon - deltaLon;
    maxLon = centerLon + deltaLon;
end

% Output confirmation to console
fprintf('Bounds: Lat[%.6f - %.6f]\n Lon[%.6f - %.6f] Spacing: %.4fm\n', minLat, maxLat, minLon, maxLon, gridSpacing);

buildingStruct = readOSM_Advanced(osmFile); % build struct so the readOSM parser can get osm infos

end


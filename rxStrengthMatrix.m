function [planeEfield, powerMatrix] = rxStrengthMatrix(rxSamplingArray, txSite, pm, gainMap)
% PLANEEFIELD calculates the Efield on each sampling point
% It uses sigstrength in chunks of rxsites to use RAM efficiently

N = numel(rxSamplingArray);

% lat = zeros(N, 1);
% lon = zeros(N, 1);

% Loop through rxSamplingArray and extract lat/lon manually
% This prevents the "Comma Separated List" explosion
% for i = 1:N
%     lat(i) = rxSamplingArray(i).Latitude;
%     lon(i) = rxSamplingArray(i).Longitude;
% end

% Try vectorized approach
lat = [rxSamplingArray.Latitude]';
lon = [rxSamplingArray.Longitude]';
numRx = numel(rxSamplingArray);

% Setup file on disk to safe progress and not lose it during a matlab crash
% :(
%fileName = 'efield_progress.mat';
fileName = 'power_progress.mat';
resultFolder = 'D:\MA_Investigation of the influence of distributed MIMO on exposure in an urban environment\MATLAB Simulation\Results';

% Build the full path
saveFile = fullfile(resultFolder, fileName);

m = matfile(saveFile, 'Writable', true);

% Initialize variable in the file (on disk, not RAM)
m.values = zeros(N, 1);

% % Pre-allocate result vector
% Let tempVals be created by sigStrength in the loop
% tempVals = zeros(N, 1);

% Define safe blocksize so the RAM doesnt overload
blockSize = 100;

% Keep track of what is happening
fprintf('Starting signal strength calculation for %d sites...\n', N);

% Initialize Output Vector
powerIncoherent_dBm = zeros(numRx, 1);

% Get Tx Power in dBm once
txPower_dBm = 10*log10(txSite.TransmitterPower) + 30;

% Isotrop Rx
rxGains = 0;

% Get Tx Mechanical Orientation
tx_mech_az = txSite.AntennaAngle(1);
tx_mech_el = txSite.AntennaAngle(2);

% Loop through blocks of rxSamplingArray
for i = 1:blockSize:N

    % Define the start and end indices for the current chunk
    idxStart = i;
    idxEnd = min(i + blockSize - 1, N); % min handles the last partial chunk
    chunkSize = idxEnd - idxStart + 1;
    % Extract the subset of receivers
    rxChunk = rxSamplingArray(idxStart:idxEnd).';

    % Ray trace chunk blockSize of receivers
    raysPerChunk = raytrace(txSite, rxChunk, pm);


    % Pre-allocated vector for results
    chunkPower_dBm = zeros(chunkSize, 1);


    % Process rays for each receiver in the current chunk
    for j = 1:chunkSize

    currentRays = raysPerChunk{j};

    % This should work already
    % Check if any rays exist (if blocked, power is -Inf)
    if isempty(currentRays)
        chunkPower_dBm(j) = -Inf;
        continue;
    end

        % Extract Data
        pathLosses  = [currentRays.PathLoss];
        AoDAngles   = [currentRays.AngleOfDeparture]; % [2 x NumRays]

        % Get azimuth and elevation angles
        AoD_az = AoDAngles(1, :);
        AoD_el = AoDAngles(2, :);

        % E. Apply Antenna Pattern
        % Convert Global Ray Angles to Local Antenna Angles
        az_local = wrapTo180(AoD_az - tx_mech_az);
        el_local = wrapTo180(AoD_el - tx_mech_el);

        % Interpolate Tx Gain using LOCAL angles
        txGains = interp2(gainMap.az, gainMap.el, gainMap.txGrid, az_local, el_local, 'linear');

        % F. Calculate Incoherent Sum
        % Net Gain per ray: Gain = G_tx + G_rx - PathLoss
        netGain_dB = txGains + rxGains - pathLosses;

        % Convert to Linear Power (Watts) relative to Tx Power (Attenuation factor)
        pathPowers_Linear = 10.^(netGain_dB ./ 10);

        % Sum the Linear Powers (Incoherent Addition)
        totalPathPower_Linear = sum(pathPowers_Linear);

        % Convert Total Linear Factor back to dB
        totalGain_dB = 10*log10(totalPathPower_Linear);

        % Final Result for this receiver
        chunkPower_dBm(j) = txPower_dBm + totalGain_dB;
    end

    % Calculate sigstrength for just this chunk
    % sigstrength does coherent addition of all rays
    % Write directly into the pre-allocated 'tempVals' vector
    % % % chunkVals = sigstrength(rxChunk, txSite, ...
    % % %     PropagationModel = pm, ...
    % % %     Type = 'power');
    % % % 
    % % % % Write small vector into the specific slot on the hard drive
    % % % % Transpose for vector dimension fit
    % % % chunkVals = chunkVals';
    % % % m.values(idxStart:idxEnd, 1) = chunkVals;

    m.values(idxStart:idxEnd, 1) = chunkPower_dBm;

    % Memory Cleanup
    % clear rxChunk chunkVals;

    clear rxChunk raysPerChunk chunkPower_dBm;

    % Print progress
    fprintf('Processed %d to %d (%.1f%%)\n', idxStart, idxEnd, (idxEnd/N)*100);

end



% % Load result back from disk to return it
storedValues = m.values;

% % Matrix of Efield values at sampling positions
planeEfield = [storedValues, lat, lon];
% planeEfield = [];

powerIncoherent_dBm = m.values;

% Create the output matrix [Power(dBm), Latitude, Longitude]
powerMatrix = [powerIncoherent_dBm, lat, lon];
% powerMatrix = [];

fprintf('Calculation complete.\n');

end
% -----------------------------------------------------------------
% Master Thesis:
% Investigation of the influence of distributed MIMO on exposure in an urban environment
% -----------------------------------------------------------------
%
% Purpose:
%   This script loads a map of the Wüllnerstraße and does a Tx-Rx DMIMO
%   simulation. The channel matrix is calculated with raytracing. A PDSCH
%   simulation is done afterwards to investigate the exposition and
%   downlink throughput.
%
%
% Author:   Christian Rieger
% Date:     04-Nov-2025
%
% Requires:
%   - 'Aachen_Wüllnerstraße.osm'
%   - Toolboxes: 5G, RF, Antenna, Communications, DSP, Mapping, Phased
%                Array System, Signal Processing
% -----------------------------------------------------------------

clc
clear
close all
tic % Time how long the simulation runs

%% Enable Distributed MIMO
% First try to implement the Massive MIMO case as reference, DMIMO as
% second case when the simulation runs, baisically add a second BS and
% change H with the information of both channels
enableDMIMO = 1;

%% Load parameter and simulation site
% When changeing map keep coordiantes from Tx-Rx in mind and change in
% parameter
fprintf('Loading Simulation Parameters, setting up Tx and Rx... \n');
osmFile = "Aachen_Wüllnerstraße.osm"; % Aachen_Wüllnerstraße.osm BauIng_Königshügel.osm AachenCity.osm
viewer = siteviewer(Buildings="Aachen_Wüllnerstraße.osm",Basemap="openstreetmap"); % Aachen_Wüllnerstraße.osm  BauIng_Königshügel.osm Weissenberg.osm


%% Load Antenna and 5G parameters
parameter;

%% Set Tx power based on allocated RBs

% getNumElements usually accounts for the total sensors, unlike prod(array.Size)
numElementsTx = getNumElements(txArray); % 128?
numElementsRx = getNumElements(rxArray); % 4?

% Calculate Gain from pattern for power per SC
weightVecTx = ones(numElementsTx, 1) / sqrt(numElementsTx);
weightVecRx = ones(numElementsRx, 1) / sqrt(numElementsRx);

% Get gain pattern for correct calculation of power in rays
[patTx, az, el] = pattern(txArray, fc, -180:180, -90:90, ...
    Weights = weightVecTx,...
    Type = 'powerdb',...
    Normalize = false);

% Get Rx pattern
[patRx, az, el] = pattern(rxArray, fc, -180:180, -90:90, ...
    Weights = weightVecRx,...
    Type = 'powerdb',...
    Normalize = false);

% [patTx, az, el] = pattern(txArray, fc, -180:180, -90:90, ...
%     Type = 'powerdb',...
%     Normalize = false);
% 
% % Get Rx pattern
% [patRx, az, el] = pattern(rxArray, fc, -180:180, -90:90, ...
%     Type = 'powerdb',...
%     Normalize = false);

gain.az = -180:180; % The X-axis (must match the pattern input)
gain.el = -90:90; % The Y-axis (must match the pattern input)
gain.txGrid = patTx; % The 3D data matrix you just calculated
gain.rxGrid = patRx; % Isotrop receiver
gain.txPeak = max(patTx, [], 'all');
gain.rxPeak = max(patRx, [], 'all'); % (dBi)

% Should EIRP also be varied on both bs?
% Ask tam

% Calculated power per subcarrier
txPowerSC = EIRP - gain.txPeak - 10 * log10(numSubcarrier); % Total EIRP - Gain - #SC (dBm)

% Calculated power per subcarrier (dBm)
% + number of allocated subcarriers (dB)
allocatedTxPower = txPowerSC + 10 * log10(numSubcarrier); % (dBm)
txPower = 10^((allocatedTxPower - 30) / 10); % Transmitted Power in W because txsite needs linear value


%% Setup transmitters and receiver
% lat, lon, azimuth
txCoords = [txLatWiwi, txLonWiwi, 90;    
            txLatAcademica, txLonAcademica, -90];

%% Set DMIMO with 1 -> both base stations otherwise just SuperC as reference
numBS = 2;

for i = 1:numBS
    txArrayCopy = clone(txArray); % For not setting the same array changes to both
    tx(i) = txsite( ...
        'Latitude', txCoords(i,1), ...
        'Longitude', txCoords(i,2), ...
        'Antenna', txArrayCopy, ...
        'AntennaAngle', [txCoords(i,3); 0], ...
        'AntennaHeight', bsHeight, ...
        'TransmitterFrequency', fc, ...
        'TransmitterPower', txPower);
end

% Specify where the first UE is located
rx = rxsite( ...
    Latitude = rxLat,...
    Longitude = rxLong,...
    AntennaHeight = ueHeight, ...
    Antenna = rxArray,...
    AntennaAngle = [0; 0]);


% Show in the siteviewer where Tx is and how the pattern looks like
show(tx);
show(rx);

if numBS > 1
    pattern(tx(1));
    pattern(tx(2));
else
    pattern(tx);
end
pattern(rx, fc, "Size", 1);

%% Setup raytracing model
% Set the propagation Model
pm = propagationModel("raytracing",... % Raytracing Parameter
    Method = "sbr",... % shooting and bouncing rays
    AngularSeparation = "low",... % "low" -> 655,362 rays launched or down to 0.05 but no memory on low Average number of degrees between launched rays, Is that good enough
    MaxRelativePathLoss = 30,... % Discard propagation paths based on a threshold [dB] TR38.901 8.4 Step 8 -25dB cutoff
    MaxNumDiffractions = 1,... % Einmal ohne probieren
    MaxNumReflections = 3);

% Calculate rays from Tx to Rx with set prop Model
fprintf('Start first raytracing calculation... \n');
rays = raytrace(tx, rx, pm, Type = "pathloss", ColorLimits = [45 180]);

% Plot all rays from both Tx
% cellfun = cell function -> OutputArray = cellfun(@FunctionName, MyCellArray);
% perform given function to each cell (like a for loop)
cellfun(@plot, rays);

% Calculate free space path loss
% L=20log10(4*pi*R/lambda)
distanceTxRx = distance(tx, rx); % (m)
freeSpacePathLoss = fspl(distanceTxRx, lambda); % (dB)

rxLOS = false(1, length(tx));
for i = 1:length(tx)
    rxLOS(i) = los(tx(i), rx);
end


% Calculate power and efield before power is focused on rx through
% beamforming
exposureBeforeBF = calcPower(tx, rx, rays, fc, allocatedTxPower);
% exposureBeforeBF2 = calcPower2(tx, rx, rays, fc, allocatedTxPower);


%% Create CDL Channel after 3GPP TR 38.901 Study on channel model for frequencies from 0.5 to 100 GHz
[channelDataBeforeBF, hEstCollectivebeforeBF] = setupDistributedChannels(rays, tx, rx, fc, UEOrientation, numRB, SCS, numSlot);
% Plot channel response in time and frequency
% figure;
% surf(pow2db(abs(hEstimation(:,:,1,1)).^2));
% shading('flat');
% xlabel('OFDM Symbols');ylabel('Subcarriers');zlabel('Magnitude Squared (dB)');
% title('Channel Magnitude Response (1^{st} tx - 1^{st} rx antenna)');


%% Beamforming/Precoding
% Using SVD
% Use numRB for channel average, 1 RB -> 12 Subcarrier, narrow band
% more RBs -> wideband simulation -> robustness over bw
[weightsTxCollective, weightsRxCollective, ~] = getBeamformingWeights(hEstCollectivebeforeBF, numLayers, scOffset, numRB);


%% Distribute weights to each BS (eq 6.3) -> slice the array and normalize per-AP power (eq 6.6)
numElementsPerBS = getNumElements(tx(1).Antenna); % Assumes same array at each BS

weightsPerBS = cell(1, numBS);

for i = 1:numBS
    startIdx = (i-1) * numElementsPerBS + 1;
    endIdx = i * numElementsPerBS;
    w = weightsTxCollective(startIdx:endIdx, :);

    % Per-BS power normalization (eq 6.6 Björnson)
    % Each BS must not exceed its own txPower budget
    w = w / norm(w, 'fro') * sqrt(numLayers);

    weightsPerBS{i} = w;
end


% Plot all beams from both transmitters
plotSitePatterns(tx, weightsPerBS, fc);

% % Do 1 layer (RANK1) or 4 layers (RANK4)
% if enableDMIMO == 1
%     tx(1).Antenna = clone(channelDataBeforeBF.BS1.Channel.TransmitAntennaArray);
%     tx(2).Antenna = clone(channelDataBeforeBF.BS2.Channel.TransmitAntennaArray);
% 
%     %%  Dein kombinierter Taper ist für Empfangsleistungs‑Patterns
%     % Use combined weights for exposure calculation
%     combinedWeightsSuperC = sum(weightsPerBS{1}(:, 1:numLayers), 2);
%     combinedWeightsAcademica = sum(weightsPerBS{2}(:, 1:numLayers), 2);
% 
%     % Jeder Singulärvektor hat bereits ||w_i||_2 = 1 und ist normiert. Aber ||combinedweights||_2
%     % ist nicht mehr normiert und wäre sqrt(2). Aber aus SVD sind w_i orthogonal daher muss nochmal
%     % normiert werden
%     combinedWeightsSuperC = combinedWeightsSuperC / norm(combinedWeightsSuperC);
%     combinedWeightsAcademica = combinedWeightsAcademica / norm(combinedWeightsAcademica);
% 
%     % Taper Tx antenna array with calculated weights (beamforming/precoding)
%     tx(1).Antenna.Taper = combinedWeightsSuperC;
%     tx(2).Antenna.Taper = combinedWeightsAcademica;
% else
%     tx(1).Antenna = clone(channelDataBeforeBF.channelSuperC.Channel.TransmitAntennaArray);
%     combinedWeightsSuperC = sum(weightsSuperC(:, 1:numLayers), 2);
%     combinedWeightsSuperC = combinedWeightsSuperC / norm(combinedWeightsSuperC);
%     tx(1).Antenna.Taper = combinedWeightsSuperC;
% end

%% Exposure: Combined taper approach (single BS and distributed precoding)
% Defensible approximation: represents worst-case instantaneous excitation
for i = 1:numBS
    tx(i).Antenna = clone(channelDataBeforeBF.("BS"+i).Channel.TransmitAntennaArray);
    combinedW = sum(weightsPerBS{i}(:, 1:numLayers), 2);
    tx(i).Antenna.Taper = combinedW / norm(combinedW);
end


% Show pattern after beamforming
if enableDMIMO == 1
    pattern(tx(1));
    pattern(tx(2));
else
    pattern(tx);
end


%% Calculate Signal strength before and after beamforming
% Do it coherent (sigstrength function) and incoherent
exposureAfterBF =  calcPower(tx, rx, rays, fc, allocatedTxPower);

% Calculate new Channel matrix after beamforming
[channelDataAfterBF, hEstCollectiveAfterBF] = setupDistributedChannels(rays, tx, rx, fc, UEOrientation, numRB, SCS, numSlot);

%% Calculate SNR and used Modulation
% Open MCS table from 38.214-Table 5.1.3.1-2
MCSIndexTablePDSCH;

% SCS should be 30kHz here -> *1e3
allocBW = 12 * SCS * 1e3 * numRB; % Bandwidth ~90MHz but its not - Guardbands? Isnt needed in calc SNR

% Calculate modulation Order, used Modulation, SNR, coderate and SE
[modOrder, usedMod, SNR_dB, coderate, SE] = calcSNRModulation(exposureAfterBF, allocBW, SCS, MCSIndexTable);


%% Run PDSCH simulation for data rate estimation
% Check if this works with LOS, NLOS and moving Rx
%
% Provide nrCDLChannel based on raytracing
% and 5G simulation parameters:
% # Frames to simulate = 10
% # MIMO layers for PDSCH = 1
% # HARQ processes = 16
% # resource blocks in the carrier grid
% coderate = Target code rate used by DL-SCH expected < 1 -> coderate/1024
% Subcarrier spacing (kHz)
% txWeights = Precoding matrix/weights applied at the transmitter across antennas
% # PRBs allocated to the PDSCH
% Starting PDSCH Simulation
% fprintf("Starting calculation of download speed: \n")
% fprintf("Start PDSCH transmission...\n")
% Theoretical throughput for 1 layer should be around 400Mbps

% Worst case -> all RB = complete bw
numRBalloc = numRB;

% for i = 1:numBS
%     release(channelDataBeforeBF.("BS"+i).Channel);
%     channelDataBeforeBF.("BS"+i).Channel.ChannelFiltering = true;
%     channelDataBeforeBF.("BS"+i).Channel.TransmitAntennaArray.Taper = 1;
%     % Mbps(i) = runPDSCH(channelDataBeforeBF.BS1.Channel, numFrames, numLayers, numHARQ, SNR_dB(1), usedMod{1}, numRB, coderate(1)/1024, SCS, weightsPerBS{1}.', numRBalloc);
%     Mbps(i) = runPDSCH2(channelDataBeforeBF.("BS"+i).Channel, numFrames, numLayers, numHARQ, SNR_dB(i), usedMod{i}, numRB, coderate(i)/1024, SCS, numRBalloc);
% end


% Prepare channels for parallel pdsch simulation sequentially
bsChannels = cell(1, numBS);
for i = 1:numBS
    ch = channelDataBeforeBF.("BS"+i).Channel;
    release(ch);
    ch.ChannelFiltering = true;
    ch.TransmitAntennaArray.Taper = 1;
    bsChannels{i} = ch;
end

% Run PDSCH simulation in parallel
Mbps = zeros(1, numBS);
parfor i = 1:numBS
    Mbps(i) = runPDSCH2(bsChannels{i}, numFrames, numLayers, numHARQ, ...
                        SNR_dB(i), usedMod{i}, numRB, coderate(i)/1024, ...
                        SCS, numRBalloc);
end


% if enableDMIMO == 1
%     % Enable Channel filter for throughput calculation
%     release(channelDataAfterBF.channelSuperC.Channel);
%     release(channelDataAfterBF.channelAcademica.Channel);
%     channelDataAfterBF.channelSuperC.Channel.ChannelFiltering = true;
%     channelDataAfterBF.channelAcademica.Channel.ChannelFiltering = true;
%     throughputMbpsSuperC = runPDSCH(channelDataAfterBF.channelSuperC.Channel, numFrames, numLayers, numHARQ, SNR_dB(1), usedMod{1}, numRB, coderate(1)/1024, SCS, weightsSuperC.', numRBalloc);
%     throughputMbpsAcademica = runPDSCH(channelDataAfterBF.channelAcademica.Channel, numFrames, numLayers, numHARQ, SNR_dB(2), usedMod{2}, numRB, coderate(1)/1024, SCS, weightsAcademica.', numRBalloc);
% 
% else
%     % release(channelDataAfterBF.channelSuperC.Channel);
%     % channelDataAfterBF.channelSuperC.Channel.ChannelFiltering = true;
%     % % channelDataAfterBF.channelSuperC.Channel.TransmitAntennaArray.Taper = 1;
%     % throughputMbpsSuperC = runPDSCH(channelDataAfterBF.channelSuperC.Channel, numFrames, numLayers, numHARQ, SNR_dB(1), usedMod{1}, numRB, coderate(1)/1024, SCS, weightsSuperC.', numRBalloc);
%     release(channelDataBeforeBF.channelSuperC.Channel);
%     channelDataBeforeBF.channelSuperC.Channel.ChannelFiltering = true;
%     % channelDataBeforeBF.channelSuperC.Channel.TransmitAntennaArray.Taper = 1;
%     throughputMbpsSuperC = runPDSCH(channelDataBeforeBF.channelSuperC.Channel, numFrames, numLayers, numHARQ, SNR_dB(1), usedMod{1}, numRB, coderate(1)/1024, SCS, weightsSuperC.', numRBalloc);
% end



%% Discretise the map into 1x1m planes of the Efield to plot a averaged Efield
% Sampling parameters set in parameter.m

% =1 Measure a 1x1m square around provided Rx coordiantes with lambda/10 sampling
% =0 Measure complete map to get an overview in 1x1m squares
measureRx = 0;
fprintf("Start calculation of Efield... \n")
% Calculate Bounds of Map (where to sample Efield)
[minLat, maxLat, minLon, maxLon, gridSpacing, buildingStruct] = getSimulationBounds(measureRx, osmFile, rxLat, rxLong, lambda);

% Calculate Grid of sampling points for each rx
[samplingPositions, gridDimension] = generateMeasurementGrid(minLat, maxLat, minLon, maxLon, gridSpacing, ueHeight, buildingStruct);

% Set a rxsite to each sampling point
% rxsite wants row vectors -> transpose
rxSamplingArray = rxsite(Latitude = samplingPositions(:, 1)',...
    Longitude = samplingPositions(:, 2)',...
    AntennaHeight = samplingPositions(:, 3)');

%% Calculate Efield at each samplepoint
% For one BS at Tx through the raytraced pathmodel pm
% [planeEfield, powerMatrix] = rxStrengthMatrix(rxSamplingArray, tx, pm, gainMap);
powerMatrix = rxStrengthMatrix2(rxSamplingArray, tx, pm, fc, allocatedTxPower);


%% Reconstruct and Plot
% Create a placeholder matrix filled with NaNs (Transparent)
% Use the dimensions we saved in gridInfo
EfieldMatrix = NaN(gridDimension.rows, gridDimension.cols);

% Fill the valid spots
% The order of 'planeEfield' matches the order of 'validMask'
% because we never re-sorted the arrays.
% planeEfield(:, 1) is the signal strength column.
EfieldMatrix(gridDimension.validMask) = planeEfield(:, 1);

%% 3D Surface Plot
figure('Color', 'w');

% SURF: X=Lon, Y=Lat, Z=SignalStrength
% This creates the "Mountain" effect where height = signal strength
hSurf = surf(gridDimension.lonGrid, gridDimension.latGrid, EfieldMatrix);
% hSurf = surf(gridDimension.lonGrid, gridDimension.latGrid, powerMatrixPlane);

set(hSurf, 'EdgeColor', 'none'); % Remove black grid lines
shading interp; % Smooth out the colors
colormap('jet'); % Rainbow colors
colorbar;
title('E-Field Strength Surface');















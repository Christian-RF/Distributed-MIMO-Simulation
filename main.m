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
% enableDMIMO = 0;
%% Set DMIMO with 1 -> both base stations otherwise just SuperC as reference
% Reference single MIMO:  numBS = 1
% Distributed MIMO:       numBS = 2
numBS = 1;

% Set up parallel workers (cores) to run PDSCH sim in parallel
% Dont use it seems to get stuck...
% pool = gcp('nocreate');
% if isempty(pool)
%     pool = parpool('local', 2);
% end
% pool.IdleTimeout = 600;  % 10 hours

outputDir = "C:\Users\rieger\Desktop\MA Christian\Results";
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

% Define resource allocations for simulation
% each will do an exposure map
rbAllocations = {
    [245, 245], "TrueDMIMO";    % Full D-MIMO
    [184, 184], "3_4";
    [123, 123], "1_2";
    [82,  82],  "1_3";
    [62,  62],  "1_4";
    [31,  31],  "1_8";};

% Fixed MCS indices to test (MATLAB 1-indexed)
fixedMCSIndices = [6, 11, 12, 20, 21, 28];

% Master results accumulator
allResults = table();
masterFile = fullfile(outputDir, "masterResults.mat");
if isfile(masterFile)
    loaded = load(masterFile, 'allResults');
    allResults = loaded.allResults;
end



%% Load parameter and simulation site
% When changeing map keep coordiantes from Tx-Rx in mind and change in
% parameter
fprintf('Loading Simulation Parameters, setting up Tx and Rx... \n');
osmFile = "Aachen_Wüllnerstraße.osm"; % Aachen_Wüllnerstraße.osm SuperC_small for testing.osm BauIng_Königshügel.osm AachenCity.osm verySmallTesting.osm
viewer = siteviewer(Buildings = osmFile, Basemap = "openstreetmap"); % Aachen_Wüllnerstraße.osm  BauIng_Königshügel.osm Weissenberg.osm


%% Load Antenna and 5G parameters
parameter;

%% Set Tx power based on allocated RBs

% getNumElements usually accounts for the total sensors, unlike prod(array.Size)
numElementsTx = getNumElements(txArray); % 64 or 128 if polarization taken into account but then calcPower breaks down but layers work better
numElementsRx = getNumElements(rxArray); % 4 or 8

% % Calculate Gain from pattern for power per SC
weightVecTx = ones(numElementsTx, 1) / sqrt(numElementsTx);
weightVecRx = ones(numElementsRx, 1) / sqrt(numElementsRx);

% [patTx, az, el] = pattern(txArray, fc, -180:180, -90:90, ...
%     Weights = weightVecTx, ...
%     Type = 'powerdb', ...
%     Normalize = false);
% Single struct, calcPower indexes gain(1) for all BSs
[patTx, az, el] = pattern(txArray, fc, -180:180, -90:90, ...
        Type = 'directivity');
    % 'Type', 'powerdb', 'Normalize', false);
gainTx.az   = -180:180;
gainTx.el   = -90:90;
gainTx.Grid = patTx;
gainTx.Peak = max(patTx, [], 'all');

% Get Rx pattern
% [patRx, ~, ~] = pattern(rxArray, fc, -180:180, -90:90, ...
%     Weights = weightVecRx, ...
%     Type = 'powerdb', ...
%     Normalize = false);
[patRx, ~, ~] = pattern(rxArray, fc, -180:180, -90:90, ...
        Type = 'directivity');
    % 'Type', 'powerdb', 'Normalize', false);
gainRx.az     = -180:180;
gainRx.el     = -90:90;
gainRx.rxGrid = patRx;
%
% gain.az = -180:180; % The X-axis (must match the pattern input)
% gain.el = -90:90; % The Y-axis (must match the pattern input)
% gain.txGrid = patTx; % The 3D data matrix you just calculated
% gain.rxGrid = patRx; % Isotrop receiver
% gain.txPeak = max(patTx, [], 'all');
% gain.rxPeak = max(patRx, [], 'all'); % (dBi)

% Should EIRP also be varied on both bs?
% Ask tam

%% RB Allocation per BS
% Define the split, change each run, give both basestation the same amount
% True DMIMO: numRB_perBS = [245, 245] == Reference: numRB_perBS = [245, 0]
% 3/4 / 3/4: numRB_perBS = [184, 184]
% 1/2 / 1/2: numRB_perBS = [123, 123]
% 1/3 / 1/3: numRB_perBS = [82, 82]
% 1/4 / 1/4: numRB_perBS = [62, 62]
% 1/8 / 1/8: numRB_perBS = [31, 31]
% numRB_perBS = [245, 245];

% Calculated power per subcarrier
txPowerSC = EIRP - gainTx.Peak - 10 * log10(numSubcarrier); % Total EIRP - Gain - #SC (dBm)

% Make it an array so calcPower can do gain(i) uniformly
gainTx = repmat(gainTx, 1, numBS);

% Calc Grid just once (stays the same) befor looping through resourceblock allocation 
% Raytracing, channel estimation, and beamforming are independent
% of RB allocation too and could be moved outside the loop for efficiency.
% Kept inside for code clarity and to avoid restructuring risk even more.
measureRx = 0;
fprintf("Start calculation of Efield... \n")
[minLat, maxLat, minLon, maxLon, gridSpacing, buildingStruct] = getSimulationBounds(measureRx, osmFile, rxLat, rxLong, lambda);
[samplingPositions, gridDimension] = generateMeasurementGrid(minLat, maxLat, minLon, maxLon, gridSpacing, ueHeight, buildingStruct);
rxSamplingArray = rxsite(Latitude = samplingPositions(:, 1)', ...
    Longitude = samplingPositions(:, 2)', ...
    AntennaHeight = samplingPositions(:, 3)');

for rbIdx = 1:size(rbAllocations, 1)
    numRB_perBS = rbAllocations{rbIdx, 1};
    rbLabel = rbAllocations{rbIdx, 2};

    % % Calculated power per subcarrier (dBm)
    % % + number of allocated subcarriers (dB)
    % allocatedTxPower = txPowerSC + 10 * log10(numSubcarrier); % (dBm)
    % txPower = 10^((allocatedTxPower - 30) / 10); % Transmitted Power in W because txsite needs linear value
    for i = 1:numBS
        if numRB_perBS(i) > 0
            numSC_i = numRB_perBS(i) * 12;
            allocatedTxPower(i) = txPowerSC + 10 * log10(numSC_i); % dBm
            txPower(i) = 10^((allocatedTxPower(i) - 30) / 10);     % Watts
        end
    end

    %% Setup transmitters and receiver
    % lat, lon, azimuth
    txCoords = [txLatWiwi, txLonWiwi, 90;
        txLatAcademica, txLonAcademica, -90];

    for i = 1:numBS
        txArrayCopy = clone(txArray); % For not setting the same array changes to both
        tx(i) = txsite( ...
            'Latitude', txCoords(i,1), ...
            'Longitude', txCoords(i,2), ...
            'Antenna', txArrayCopy, ...
            'AntennaAngle', [txCoords(i,3); 0], ...
            'AntennaHeight', bsHeight, ...
            'TransmitterFrequency', fc, ...
            'TransmitterPower', txPower(i));
    end

    % Specify where the first UE is located
    rx = rxsite( ...
        Latitude = rxLat,...
        Longitude = rxLong,...
        AntennaHeight = ueHeight, ...
        Antenna = rxArray,...
        AntennaAngle = [0; 0]);


    % Show in the siteviewer where Tx is and how the pattern looks like
    % % show(tx);
    % % show(rx);
    % %
    % % if numBS > 1
    % %     pattern(tx(1));
    % %     pattern(tx(2));
    % % else
    % %     pattern(tx);
    % % end
    % % pattern(rx, fc, "Size", 1);

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
    % cellfun(@plot, rays);

    % Calculate free space path loss
    % L=20log10(4*pi*R/lambda)
    distanceTxRx = distance(tx, rx); % (m)
    freeSpacePathLoss = fspl(distanceTxRx, lambda); % (dB)

    rxLOS = false(1, length(tx));
    for i = 1:numBS
        rxLOS(i) = los(tx(i), rx);
    end


    % Calculate power and efield before power is focused on rx through
    % beamforming
    % exposureBeforeBF = calcPower(tx, rx, rays, fc, allocatedTxPower);
    % exposure = calcPower3(tx, rx, rays, fc, allocatedTxPower, gainTx, gainRx);
    % exposure_sigstrength = sigstrength(rx, tx, pm, Type = 'power')
    % exposure_pattern =  calcPower(tx, rx, rays, fc, allocatedTxPower)
    % exposure_interp2 = calcPower3(tx, rx, rays, fc, allocatedTxPower, gainTx, gainRx)
    exposure = calcPower4(tx, rx, rays, fc, allocatedTxPower, gainTx, gainRx);
    % for i = 1:numBS
    %     exposure(i) = calcPower3(tx(i), rx, rays(i), fc, allocatedTxPower(i), gainTx(i), gainRx);
    % end


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

    weightsPerBS = cell(1, numBS);

    for i = 1:numBS
        startIdx = (i-1) * numElementsTx + 1;
        endIdx = i * numElementsTx;
        w = weightsTxCollective(startIdx:endIdx, :);

        % Per-BS power normalization (eq 6.6 Björnson)
        % Each BS must not exceed its own txPower budget
        w = w / norm(w, 'fro') * sqrt(numLayers);

        weightsPerBS{i} = w;
    end


    % % % % % Plot all beams from both transmitters
    % % % % plotSitePatterns(tx, weightsPerBS, fc);

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
    % Defensible approximation: represents worst-case instantaneous
    % exposure
    for i = 1:numBS
        tx(i).Antenna = clone(channelDataBeforeBF.("BS"+i).Channel.TransmitAntennaArray);
        combinedW = sum(weightsPerBS{i}(:, 1:numLayers), 2);
        tx(i).Antenna.Taper = combinedW / norm(combinedW);

        % For faster gain calculation
        [patTxBF_i, az, el] = pattern(tx(i).Antenna, fc, -180:180, -90:90, ...
            Type = 'directivity');
            % 'Type', 'powerdb', 'Normalize', false);
        gainTxBF(i).az = -180:180;
        gainTxBF(i).el = -90:90;
        gainTxBF(i).Grid = patTxBF_i;
        gainTxBF(i).Peak = max(patTxBF_i, [], 'all');
    end

    % % % Show pattern after beamforming
    % % if numBS == 2
    % %     pattern(tx(1));
    % %     pattern(tx(2));
    % % else
    % %     pattern(tx);
    % % end

    % dynRange = 40; % dynamic range to plot pattern similar to pattern()
    % figure('Color', 'w');
    % for i = 1:numBS
    %     [AZ, EL] = meshgrid(deg2rad(gainTxBF(i).az), deg2rad(gainTxBF(i).el));
    %     pat = gainTxBF(i).Grid;
    % 
    %     patClipped = max(pat, gainTxBF(i).Peak - dynRange);
    %     R = (patClipped - (gainTxBF(i).Peak - dynRange)) / dynRange;
    % 
    %     X = R .* cos(EL) .* cos(AZ);
    %     Y = R .* cos(EL) .* sin(AZ);
    %     Z = R .* sin(EL);
    % 
    %     subplot(1, numBS, i);
    %     surf(X, Y, Z, patClipped, 'EdgeColor', 'none');
    %     colormap(jet);
    %     caxis([gainTxBF(i).Peak - dynRange, gainTxBF(i).Peak]);
    %     colorbar;
    %     axis equal;
    %     title(sprintf('BS %d (%.1f dBi)', i, gainTxBF(i).Peak));
    % end


    %% Calculate Signal strength before and after beamforming
    % Do it coherent (sigstrength function) and incoherent
    % exposureAfterBF =  calcPower3(tx, rx, rays, fc, allocatedTxPower, gainTxBF, gainRx);

    % exposureAfterBF_sigstrength = sigstrength(rx, tx, pm, Type = 'power');
    % exposureAfterBF_pattern = calcPower(tx, rx, rays, fc, allocatedTxPower);
    % exposureAfterBF_interp2 = calcPower3(tx, rx, rays, fc, allocatedTxPower, gainTxBF, gainRx);
    exposureAfterBF = calcPower4(tx, rx, rays, fc, allocatedTxPower, gainTxBF, gainRx);
    % for i = 1:numBS
    %     exposureAfterBF(i) =  calcPower3(tx(i), rx, rays(i), fc, allocatedTxPower(i), gainTxBF(i), gainRx);
    % end

    %% Discretise the map into 1x1m planes of the Efield to plot a averaged Efield
    % Sampling parameters set in parameter.m

    % =1 Measure a 1x1m square around provided Rx coordiantes with lambda/10 sampling
    % =0 Measure complete map to get an overview in 1x1m squares
    % measureRx = 0;
    % fprintf("Start calculation of Efield... \n")
    % % Calculate Bounds of Map (where to sample Efield)
    % [minLat, maxLat, minLon, maxLon, gridSpacing, buildingStruct] = getSimulationBounds(measureRx, osmFile, rxLat, rxLong, lambda);
    % 
    % % Calculate Grid of sampling points for each rx
    % [samplingPositions, gridDimension] = generateMeasurementGrid(minLat, maxLat, minLon, maxLon, gridSpacing, ueHeight, buildingStruct);
    % 
    % % Set a rxsite to each sampling point
    % % rxsite wants row vectors -> transpose
    % rxSamplingArray = rxsite(Latitude = samplingPositions(:, 1)',...
    %     Longitude = samplingPositions(:, 2)',...
    %     AntennaHeight = samplingPositions(:, 3)');

    %% Calculate Efield at each samplepoint
    % For one BS at Tx through the raytraced pathmodel pm
    % [planeEfield, powerMatrix] = rxStrengthMatrix(rxSamplingArray, tx, pm, gainMap);
    % powerMatrix = rxStrengthMatrix2(viewer, rxSamplingArray, tx, pm, fc, allocatedTxPower);
    % powerMatrix = rxStrengthMatrix3(viewer, rxSamplingArray, tx, pm, fc, allocatedTxPower, gainTx, gainRx);
    %% Try to just calc powerMatrix once
    %%%%% powerMatrix = rxStrengthMatrix3(viewer, rxSamplingArray, tx, pm, fc, allocatedTxPower, gainTxBF, gainRx);


    % Compute exposure map only on the FIRST RB iteration
    % Since it stays the same, same pattern, same map -> channel, same
    % position just power changes via resource allocation
    % Calc constant offset from map
    if rbIdx == 1
        %% Calc powerMatrix
        powerMatrix_ref = rxStrengthMatrix4(viewer, rxSamplingArray, tx, pm, fc, allocatedTxPower, gainTxBF, gainRx);
        allocatedTxPower_ref = allocatedTxPower;
        save(fullfile(outputDir, sprintf("MP%d_referenceMap.mat", mpIdx)), ...
            'powerMatrix_ref', 'allocatedTxPower_ref', 'gridDimension');
        powerMatrix = powerMatrix_ref;
    else
        powerMatrix = powerMatrix_ref;

        % Per-BS scaling
        for i = 1:numBS
            delta_dB = allocatedTxPower(i) - allocatedTxPower_ref(i);

            % Power fields: + delta_dB
            powerMatrix.coherentPower_dBm(:, i)      = powerMatrix_ref.coherentPower_dBm(:, i) + delta_dB;
            powerMatrix.incoherentPower_dBm(:, i)    = powerMatrix_ref.incoherentPower_dBm(:, i) + delta_dB;
            powerMatrix.coherentPowerIso_dBm(:, i)   = powerMatrix_ref.coherentPowerIso_dBm(:, i) + delta_dB;
            powerMatrix.incoherentPowerIso_dBm(:, i) = powerMatrix_ref.incoherentPowerIso_dBm(:, i) + delta_dB;

            % E-field fields: + delta_dB / 2 (E scales as sqrt of power)
            powerMatrix.coherentEfield_dBuv(:, i)    = powerMatrix_ref.coherentEfield_dBuv(:, i) + delta_dB / 2;
            powerMatrix.incoherentEfield_dBuv(:, i)  = powerMatrix_ref.incoherentEfield_dBuv(:, i) + delta_dB / 2;
        end

        % Totals: recompute from scaled per-BS values (not just a simple offset)
        % Incoherent total: sum powers in linear
        incPow_lin = 10.^(powerMatrix.incoherentPower_dBm / 10);  % N x numBS
        powerMatrix.totalIncoherentPower_dBm = 10*log10(sum(incPow_lin, 2));

        cohPow_lin = 10.^(powerMatrix.coherentPower_dBm / 10);
        powerMatrix.totalCoherentPower_dBm = 10*log10(sum(cohPow_lin, 2));

        % Incoherent E-field total: sum powers (isotropic), then convert
        incPowIso_lin = 10.^(powerMatrix.incoherentPowerIso_dBm / 10);
        totalIncPowIso_dBm = 10*log10(sum(incPowIso_lin, 2));
        Z0_term = 10*log10(rfprop.Constants.Z0 * 4 * pi);
        lambda_term = 20*log10(rfprop.Constants.LightSpeed / fc);
        powerMatrix.totalIncoherentEfield_dBuv = totalIncPowIso_dBm - 30 + Z0_term + 120 - lambda_term;

        cohPowIso_lin = 10.^(powerMatrix.coherentPowerIso_dBm / 10);
        totalCohPowIso_dBm = 10*log10(sum(cohPowIso_lin, 2));
        powerMatrix.totalCoherentEfield_dBuv = totalCohPowIso_dBm - 30 + Z0_term + 120 - lambda_term;
    end


    powerMatrixPlane = NaN(gridDimension.rows, gridDimension.cols);
    powerMatrixPlane(gridDimension.validMask) = powerMatrix.totalIncoherentPower_dBm;

    efieldTest = NaN(gridDimension.rows, gridDimension.cols);
    efieldTest(gridDimension.validMask) = powerMatrix.totalCoherentEfield_dBuv;
    % %% 3D Surface Plot
    fig = figure('Color', 'w');
    %
    % SURF: X=Lon, Y=Lat, Z=SignalStrength
    % This creates the "Mountain" effect where height = signal strength
    % hSurf = surf(gridDimension.lonGrid, gridDimension.latGrid, EfieldMatrix);
    hSurf = surf(gridDimension.lonGrid, gridDimension.latGrid, powerMatrixPlane);
    % hSurf = surf(gridDimension.lonGrid, gridDimension.latGrid, efieldTest);

    set(hSurf, 'EdgeColor', 'none'); % Remove black grid lines
    shading interp; % Smooth out the colors
    colormap('jet'); % Rainbow colors
    colorbar;
    savefig(fig, fullfile(outputDir, sprintf("MP%d_%s_exposureMap.fig", mpIdx, rbLabel)));
    close(fig);

    % Calculate new Channel matrix after beamforming
    % [channelDataAfterBF, hEstCollectiveAfterBF] = setupDistributedChannels(rays, tx, rx, fc, UEOrientation, numRB, SCS, numSlot);

    %% Calculate SNR and used Modulation
    % Open MCS table from 38.214-Table 5.1.3.1-2
    MCSIndexTablePDSCH;
    % SCS should be 30kHz here -> *1e3
    % allocBW = 12 * SCS * 1e3 * numRB; % Bandwidth ~90MHz but its not - Guardbands? Isnt needed in calc SNR

    % Calculate modulation Order, used Modulation, SNR, coderate and SE
    % [modOrder, usedMod, SNR_dB, coderate, SE] = calcSNRModulation(exposureAfterBF, allocBW, SCS, MCSIndexTable);
    % allocBW = 12 * SCS * 1e3 * numRB_perBS;
    % [modOrder, usedMod, SNR_dB, coderate, SE] = calcSNRModulation(exposureAfterBF_interp2, allocBW, SCS, MCSIndexTable);

    modOrder = zeros(1, numBS);
    usedMod = cell(1, numBS);
    SNR_Lin = zeros(1, numBS);
    SNR_dB = zeros(1, numBS);
    coderate = zeros(1, numBS);
    SE = zeros(1, numBS);
    mcsIdx = zeros(1, numBS);
    for i = 1:numBS
        allocBW_i = 12 * SCS * 1e3 * numRB_perBS(i);
        exposure_i.incoherentPower_dBm = exposureAfterBF.incoherentPower_dBm(i);
        [modOrder(i), usedMod{i}, SNR_dB(i), coderate(i), SE(i), mcsIdx(i)] = ...
            calcSNRModulation(exposure_i, allocBW_i, SCS, MCSIndexTable);
        SNR_Lin(i) = 10^(SNR_dB(i) / 10);
    end
    SNRjoint_dB = 10*log10(sum(SNR_Lin));

    % Prepare channels for parallel pdsch simulation sequentially
    bsChannels = cell(1, numBS);
    for i = 1:numBS
        ch = channelDataBeforeBF.("BS"+i).Channel;
        release(ch);
        ch.ChannelFiltering = true;
        ch.TransmitAntennaArray.Taper = 1;
        bsChannels{i} = ch;
    end

    % Compute PRB allocation per BS
    % Compute PRB allocation per BS
    prbSets = cell(1, numBS);
    for i = 1:numBS
        if numRB_perBS(i) > 0
            prbSets{i} = 0:(numRB_perBS(i) - 1);  % always start at 0
        else
            prbSets{i} = [];
        end
    end
    % SNR for DMIMO?
    % SNR1_lin = 10^(SNR_dB(1) / 10);
    % SNR2_lin = 10^(SNR_dB(2) / 10);
    % SNRjoint_dB = 10*log10(SNR_lin(1) + SNR_lin(2));

    % Store selected MCS indices based on SNR
    autoMCSIdx = mcsIdx;

    % Build list: auto + 6 fixed = 7 variations to loop
    mcsVariations = [NaN, fixedMCSIndices];  % NaN = marker for "auto"

    for mcsVarIdx = 1:length(mcsVariations)

        % MCS lookup
        if isnan(mcsVariations(mcsVarIdx))
            % Use auto-selected MCS (already computed above)
            currentMCS = autoMCSIdx;
            mcsLabel = "auto";
        else
            % Override all BS to this fixed MCS index
            currentMCS = repmat(mcsVariations(mcsVarIdx), 1, numBS);
            mcsLabel = sprintf("fixed%d", mcsVariations(mcsVarIdx));
        end

        % Get mod/cr from MCS table
        for i = 1:numBS
            mcsIdx(i) = currentMCS(i);
            modOrder(i) = MCSIndexTable.Modulation_Order(mcsIdx(i));
            usedMod{i} = string(MCSIndexTable.Modulation(mcsIdx(i)));
            coderate(i) = MCSIndexTable.Target_Code_Rate(mcsIdx(i));
            SE(i) = MCSIndexTable.Spectral_efficiency(mcsIdx(i)); % chose from set MCS not based on SNR
            % SNR stays the same —> depends on power, not MCS
        end

        %% Run PDSCH simulation for data rate estimation
        % Check if this works with LOS, NLOS and moving Rx
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
        % numRBalloc = numRB;

        % for i = 1:numBS
        %     release(channelDataBeforeBF.("BS"+i).Channel);
        %     channelDataBeforeBF.("BS"+i).Channel.ChannelFiltering = true;
        %     channelDataBeforeBF.("BS"+i).Channel.TransmitAntennaArray.Taper = 1;
        %     % Mbps(i) = runPDSCH(channelDataBeforeBF.BS1.Channel, numFrames, numLayers, numHARQ, SNR_dB(1), usedMod{1}, numRB, coderate(1)/1024, SCS, weightsPerBS{1}.', numRBalloc);
        %     Mbps(i) = runPDSCH2(channelDataBeforeBF.("BS"+i).Channel, numFrames, numLayers, numHARQ, SNR_dB(i), usedMod{i}, numRB, coderate(i)/1024, SCS, numRBalloc);
        % end




        % Run PDSCH simulation
        % Simulating Dual Connectivity (incoherent) rather than Joint Transmission (DMIMO) (coherent)
        % dual connectivity would split the frequency band and serve the UE but in
        % true DMIMO the signal is the same from both Tx and would add coherently
        % at UE -> +3dB power in the best case for two Tx, check how much that
        % changes the MCS probably not even worth it for two far away BS
        % For true implementation Htotal​=[H1​,H2​] and one run of PDSCH but no
        % time/knowledge to implement this with two tx into the matlab functions
        % for 5G
        Mbps_perBS = zeros(1, numBS);
        % Is parfor a good idea seems to get stuck -> not worth it
        for i = 1:numBS
            if ~isempty(prbSets{i})
                Mbps_perBS(i) = runPDSCH2(bsChannels{i}, numFrames, numLayers, numHARQ, ...
                    SNR_dB(i), usedMod{i}, numRB, coderate(i)/1024, ...
                    SCS, prbSets{i});  % pass PRB set instead of numRBalloc
            end
        end
        Mbps_total = sum(Mbps_perBS);

        %% Build result rows
        for i = 1:numBS
            row = table( ...
                mpIdx, rxLat, rxLong, i, numBS, ...
                numRB_perBS(i), string(rbLabel), ...
                string(mcsLabel), mcsIdx(i), string(usedMod{i}), coderate(i), ...
                SNR_dB(i), SNRjoint_dB, SE(i), ...
                allocatedTxPower(i), txPower(i), ...
                rxLOS(i), distanceTxRx(i), freeSpacePathLoss(i), ...
                gainTxBF(i).Peak, ...
                ... % Per-BS pre-BF
                exposure.incoherentPower_dBm(i), ...
                exposure.coherentPower_dBm(i), ...
                exposure.incoherentPowerIso_dBm(i), ...
                exposure.coherentPowerIso_dBm(i), ...
                exposure.incoherentEfield_dBuv(i), ...
                exposure.coherentEfield_dBuv(i), ...
                ... % Total pre-BF
                exposure.totalIncoherentPower_dBm, ...
                exposure.totalCoherentPower_dBm, ...
                exposure.totalIncoherentEfield_dBuv, ...
                exposure.totalCoherentEfield_dBuv, ...
                ... % Per-BS post-BF
                exposureAfterBF.incoherentPower_dBm(i), ...
                exposureAfterBF.coherentPower_dBm(i), ...
                exposureAfterBF.incoherentPowerIso_dBm(i), ...
                exposureAfterBF.coherentPowerIso_dBm(i), ...
                exposureAfterBF.incoherentEfield_dBuv(i), ...
                exposureAfterBF.coherentEfield_dBuv(i), ...
                ... % Total post-BF
                exposureAfterBF.totalIncoherentPower_dBm, ...
                exposureAfterBF.totalCoherentPower_dBm, ...
                exposureAfterBF.totalIncoherentEfield_dBuv, ...
                exposureAfterBF.totalCoherentEfield_dBuv, ...
                ... % Throughput
                Mbps_perBS(i), sum(Mbps_perBS), ...
                'VariableNames', { ...
                'mpIdx','rxLat','rxLong','BS','numBS', ...
                'numRB','rbLabel', ...
                'mcsLabel','mcsIdx','modulation','coderate', ...
                'SNR_dB','SNRjoint_dB','SE', ...
                'txPower_dBm','txPower_W', ...
                'LOS','distance_m','FSPL_dB', ...
                'BF_peakDirectivity_dBi', ...
                ... % Per-BS pre-BF
                'perBS_preBF_incohPow_dBm', ...
                'perBS_preBF_cohPow_dBm', ...
                'perBS_preBF_incohPowIso_dBm', ...
                'perBS_preBF_cohPowIso_dBm', ...
                'perBS_preBF_incohEf_dBuv', ...
                'perBS_preBF_cohEf_dBuv', ...
                ... % Total pre-BF
                'total_preBF_incohPow_dBm', ...
                'total_preBF_cohPow_dBm', ...
                'total_preBF_incohEf_dBuv', ...
                'total_preBF_cohEf_dBuv', ...
                ... % Per-BS post-BF
                'perBS_postBF_incohPow_dBm', ...
                'perBS_postBF_cohPow_dBm', ...
                'perBS_postBF_incohPowIso_dBm', ...
                'perBS_postBF_cohPowIso_dBm', ...
                'perBS_postBF_incohEf_dBuv', ...
                'perBS_postBF_cohEf_dBuv', ...
                ... % Total post-BF
                'total_postBF_incohPow_dBm', ...
                'total_postBF_cohPow_dBm', ...
                'total_postBF_incohEf_dBuv', ...
                'total_postBF_cohEf_dBuv', ...
                ... % Throughput
                'Mbps_perBS','Mbps_total'});
            allResults = [allResults; row];
        end

    end

    %% Save matrices per RB allocation (not per MCS)
    matrixFile = fullfile(outputDir, sprintf("MP%d_%s_matrices.mat", mpIdx, rbLabel));
    save(matrixFile, 'powerMatrix', 'powerMatrixPlane', 'exposure', 'exposureAfterBF');

end
close(viewer);

% Set data name depending on settings
timestamp = char(datetime, 'ddMMyyyy_HHmm');
simParameters = sprintf("MP%d_", mpIdx);
for i = 1:numBS
    simParameters = simParameters + sprintf("BS%d_RB%d_MCS%d_%s_CR%d", ...
        i, numRB_perBS(i), mcsIdx(i), string(usedMod{i}), coderate(i));
    if i < numBS
        simParameters = simParameters + "_";
    end
end

%% Save master table
save(masterFile, 'allResults');
writetable(allResults, fullfile(outputDir, "masterResults.csv"));

%% Save full workspace
workspaceFile = fullfile(outputDir, simParameters + "_" + timestamp + "_workspace.mat");
save(workspaceFile, '-v7.3');

% Save Matrices
% matrixFile = fullfile(outputDir, simParameters + "_" + timestamp + "_matrices.mat");
% save(matrixFile, 'powerMatrix', 'powerMatrixPlane', 'exposure', 'exposureAfterBF');



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




% %% Reconstruct and Plot
% % Create a placeholder matrix filled with NaNs (Transparent)
% % Use the dimensions we saved in gridInfo
% EfieldMatrix = NaN(gridDimension.rows, gridDimension.cols);

% % Fill the valid spots
% % The order of 'planeEfield' matches the order of 'validMask'
% % because we never re-sorted the arrays.
% % planeEfield(:, 1) is the signal strength column.
% EfieldMatrix(gridDimension.validMask) = planeEfield(:, 1);
%
% %% 3D Surface Plot
% figure('Color', 'w');
% %
% % SURF: X=Lon, Y=Lat, Z=SignalStrength
% % This creates the "Mountain" effect where height = signal strength
% % hSurf = surf(gridDimension.lonGrid, gridDimension.latGrid, EfieldMatrix);
% hSurf = surf(gridDimension.lonGrid, gridDimension.latGrid, powerMatrixPlane);
% % hSurf = surf(gridDimension.lonGrid, gridDimension.latGrid, efieldTest);
%
% set(hSurf, 'EdgeColor', 'none'); % Remove black grid lines
% shading interp; % Smooth out the colors
% colormap('jet'); % Rainbow colors
% colorbar;


% clim([-140 -40]);
% title('E-Field Strength Surface');


% figure('Color', 'w');
% hold on;
% % Building footprints
% for i = 1:length(buildingStruct)
%     fill(buildingStruct(i).Lon, buildingStruct(i).Lat, [0.6 0.6 0.6], ...
%          'EdgeColor', 'k');
% end
% % Grid points
% latValid = gridDimension.latGrid(gridDimension.validMask);
% lonValid = gridDimension.lonGrid(gridDimension.validMask);
% scatter(lonValid, latValid, 1, 'b', 'filled');
% axis equal; grid on;
% xlabel('Longitude'); ylabel('Latitude');
% title(sprintf('Measurement Grid (%.1fm spacing, %d points)', gridSpacing, sum(gridDimension.validMask(:))));


% figure('Color', 'w');
% hSurf = surf(gridDimension.lonGrid, gridDimension.latGrid, powerMatrixPlane);
% set(hSurf, 'EdgeColor', 'none');
% shading interp; colormap('jet');
% hold on;
%
% % Overlay actual building polygons — always smooth
% for i = 1:length(buildingStruct)
%     fill3(buildingStruct(i).Lon, buildingStruct(i).Lat, ...
%         max(powerMatrixPlane(:)) * ones(size(buildingStruct(i).Lat)), ...
%         [0.4 0.4 0.4], 'EdgeColor', 'k', 'FaceAlpha', 0.9);
% end



% %% 1) Save complete workspace (as requested by supervisor)
% workspaceFile = fullfile(outputDir, simParameters + "_" + timestamp + "_workspace.mat");
% save(workspaceFile, '-v7.3'); % -v7.3 support for large files
%
% %% 2) Build flat results table (one row per BS, appendable across runs)
% % Only scalar/short values — no large matrices here
% resultRows = cell(numBS, 1);
% for i = 1:numBS
%     resultRows{i} = table( ...
%         string(simParameters), mpIdx, rxLat, rxLong, i, numBS, ...
%         numRB_perBS(i), string(rbLabel), mcsIdx(i), string(usedMod{i}), coderate(i), ...
%         SNR_dB(i), SNRjoint_dB, SE(i), allocatedTxPower(i), txPower(i), ...
%         rxLOS(i), distanceTxRx(i), freeSpacePathLoss(i), ...
%         gainTxBF(i).Peak, ...
%         exposure.totalIncoherentPower_dBm, ...
%         exposure.totalCoherentPower_dBm, ...
%         exposureAfterBF.totalIncoherentPower_dBm, ...
%         exposureAfterBF.totalCoherentPower_dBm, ...
%         Mbps_perBS(i), sum(Mbps_perBS), ...
%         'VariableNames', { ...
%         'SimulationRun', 'mpIdx','rxLat','rxLong','BS','numBS', ...
%         'numRB','rbLabel','mcsIdx','modulation','coderate', ...
%         'SNR_dB','SNRjoint_dB','SE','txPower_dBm','txPower_W', ...
%         'LOS','distance_m','FSPL_dB', ...
%         'BF_peakDirectivity_dBi', ...
%         'totalExp_preBF_incoh_dBm', ...
%         'totalExp_preBF_coh_dBm', ...
%         'totalExp_postBF_incoh_dBm', ...
%         'totalExp_postBF_coh_dBm', ...
%         'Mbps_perBS','Mbps_total'});
% end
% resultsTable = vertcat(resultRows{:});
%
% % Append to master results
% masterFile = fullfile(outputDir, "masterResults.mat");
% if isfile(masterFile)
%     loaded = load(masterFile, 'allResults');
%     allResults = [loaded.allResults; resultsTable];
% else
%     allResults = resultsTable;
% end
% save(masterFile, 'allResults');
%

% % % Analysis of data
% % load('summaryTable.mat');
% %
% % % Compare throughput across RB splits for MP1
% % mp1 = resultsTable(resultsTable.mpIdx == 1, :);
% % figure; bar(categorical(mp1.rbLabel), mp1.Mbps_total);
% %
% % % Exposure vs throughput scatter
% % figure; scatter(resultsTable.postBF_incohPwr, resultsTable.Mbps_total);












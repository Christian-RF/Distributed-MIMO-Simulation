function [powerMatrix] = rxStrengthMatrix3(viewer, rxSamplingArray, txSite, pm, fc, allocatedTxPower, gainTx, gainRx)

N = numel(rxSamplingArray);
lat = [rxSamplingArray.Latitude]';
lon = [rxSamplingArray.Longitude]';

% Setup incremental save
% D:\MA_Investigation of the influence of distributed MIMO on exposure in an urban environment\MATLAB Simulation\Results
% C:\Users\chris\OneDrive - Students RWTH Aachen University\RWTH\Master Thesis\Results
resultFolder = 'D:\MA_Investigation of the influence of distributed MIMO on exposure in an urban environment\MATLAB Simulation\Results';
saveFile = fullfile(resultFolder, 'power_progress_1BS_fullRB_layer1_MP1.mat');
m = matfile(saveFile, 'Writable', true);
m.incoherentPower_dBm = -inf(N, 1);
m.coherentPower_dBm = -inf(N, 1);
m.incoherentEfield_dBuv = -inf(N, 1);
m.coherentEfield_dBuv = -inf(N, 1);

% Blocks of rx in sampling grid so the RAM/siteviewer crashes arent a problem 
blockSize = 500;

fprintf('Starting power calculation for %d points...\n', N);
totalBlocks = ceil(N / blockSize);
tStart = tic;

for block = 1:blockSize:N
    tBlock = tic; 
    
    idxStart = block;
    idxEnd   = min(block + blockSize - 1, N);
    chunkSize = idxEnd - idxStart + 1;
    blockNum = ceil(block / blockSize);

    % Check viewer health because ray tracer just runs properly with
    % siteviewer open??????
    if ~isvalid(viewer)
        warning('Siteviewer crashed at block %d, reopening...', block);
        viewer = siteviewer("Buildings", "Aachen_Wüllnerstraße.osm", ...
            "Basemap", "openstreetmap");
        pause(2);
    end

    rxChunk = rxSamplingArray(idxStart:idxEnd);

    % Ray trace this is the bottleneck 
    raysAllTxRx = raytrace(txSite, rxChunk, pm, Type="pathloss", Map=viewer); 

    % Diagnostic on first block
    if block == 1
        nMultipath = 0; nTotal = 0;
        for dbg = 1:numel(raysAllTxRx)
            if ~isempty(raysAllTxRx{dbg})
                for r = 1:numel(raysAllTxRx{dbg})
                    nTotal = nTotal + 1;
                    if raysAllTxRx{dbg}(r).NumReflections > 0 || ...
                            raysAllTxRx{dbg}(r).NumDiffractions > 0
                        nMultipath = nMultipath + 1;
                    end
                end
            end
        end
        fprintf('Block 1: %d/%d rays have multipath\n', nMultipath, nTotal);
        if nMultipath == 0
            warning('NO multipath — buildings may not be loaded!');
        end
    end

    incPow = -inf(chunkSize, 1);
    cohPow = -inf(chunkSize, 1);
    incEf = -inf(chunkSize, 1);
    cohEf = -inf(chunkSize, 1);

    % Simple for loop — interp2 is fast, no need for parfor
    for j = 1:chunkSize
        raysForRx = raysAllTxRx(:, j);
        hasRays = ~cellfun(@isempty, raysForRx);

        if ~any(hasRays)
            continue;
        end

        raysValid = raysForRx(hasRays);
        txValid   = txSite(hasRays);
        gainValid = gainTx(hasRays);  % matching gain grids

        if isscalar(allocatedTxPower)
            allocPowerValid = allocatedTxPower;
        else
            allocPowerValid = allocatedTxPower(hasRays);
        end

        exp = calcPower3(txValid, rxChunk(j), raysValid, fc, allocPowerValid, gainValid, gainRx);

        incPow(j) = exp.totalIncoherentPower_dBm;
        cohPow(j) = exp.totalCoherentPower_dBm;
        incEf(j) = exp.totalIncoherentEfield_dBuv;
        cohEf(j) = exp.totalCoherentEfield_dBuv;
    end

    % Save block to disk
    m.incoherentPower_dBm(idxStart:idxEnd, 1) = incPow;
    m.coherentPower_dBm(idxStart:idxEnd, 1) = cohPow;
    m.incoherentEfield_dBuv(idxStart:idxEnd, 1) = incEf;
    m.coherentEfield_dBuv(idxStart:idxEnd, 1) = cohEf;

    % Timing and progress
    blockTime = toc(tBlock);
    elapsed = toc(tStart);
    pctDone = idxEnd / N * 100;
    avgPerPoint = elapsed / idxEnd;
    remaining = avgPerPoint * (N - idxEnd);

    fprintf('Block %d/%d | %d–%d of %d (%.1f%%) | Block: %.1fs | Elapsed: %.1fm | ETA: %.1fm\n', ...
             blockNum, totalBlocks, idxStart, idxEnd, N, pctDone, ...
             blockTime, elapsed/60, remaining/60);
end

% Build output
powerMatrix.incoherentPower_dBm = m.incoherentPower_dBm;
powerMatrix.coherentPower_dBm = m.coherentPower_dBm;
powerMatrix.incoherentEfield_dBuv = m.incoherentEfield_dBuv;
powerMatrix.coherentEfield_dBuv = m.coherentEfield_dBuv;
powerMatrix.lat = lat;
powerMatrix.lon = lon;

fprintf('Power calculation complete.\n');
end
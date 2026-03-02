function [powerMatrix] = rxStrengthMatrix4(viewer, rxSamplingArray, txSite, pm, fc, allocatedTxPower, gainTx, gainRx)

N = numel(rxSamplingArray);
numBS = numel(txSite);
lat = [rxSamplingArray.Latitude]';
lon = [rxSamplingArray.Longitude]';

resultFolder = 'C:\Users\rieger\Desktop\MA Christian\Results';
saveFile = fullfile(resultFolder, 'power_progress.mat');
m = matfile(saveFile, 'Writable', true);

% Totals (4 values per point — same as before)
m.totalIncoherentPower_dBm = -inf(N, 1);
m.totalCoherentPower_dBm   = -inf(N, 1);
m.totalIncoherentEfield_dBuv = -inf(N, 1);
m.totalCoherentEfield_dBuv   = -inf(N, 1);

% Per-BS (6 values per BS per point)
m.coherentPower_dBm      = -inf(N, numBS);
m.incoherentPower_dBm    = -inf(N, numBS);
m.coherentPowerIso_dBm   = -inf(N, numBS);
m.incoherentPowerIso_dBm = -inf(N, numBS);
m.coherentEfield_dBuv    = -inf(N, numBS);
m.incoherentEfield_dBuv  = -inf(N, numBS);

m.phasorSum_perBS    = complex(zeros(N, numBS)); % complex, so zeros not -inf
m.phasorSumIso_perBS = complex(zeros(N, numBS));


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

    if ~isvalid(viewer)
        warning('Siteviewer crashed at block %d, reopening...', block);
        viewer = siteviewer("Buildings", "Aachen_Wüllnerstraße.osm", ...
            "Basemap", "openstreetmap");
        pause(2);
    end

    rxChunk = rxSamplingArray(idxStart:idxEnd);
    raysAllTxRx = raytrace(txSite, rxChunk, pm, Type="pathloss", Map=viewer);

    % Preallocate Chunk arrays 
    tot_incPow = -inf(chunkSize, 1);
    tot_cohPow = -inf(chunkSize, 1);
    tot_incEf  = -inf(chunkSize, 1);
    tot_cohEf  = -inf(chunkSize, 1);

    % PreallocateChunk arrays per BS
    bs_cohPow    = -inf(chunkSize, numBS);
    bs_incPow    = -inf(chunkSize, numBS);
    bs_cohPowIso = -inf(chunkSize, numBS);
    bs_incPowIso = -inf(chunkSize, numBS);
    bs_cohEf     = -inf(chunkSize, numBS);
    bs_incEf     = -inf(chunkSize, numBS);

    % Preallocate phasors for scaling
    bs_phasorSum    = complex(zeros(chunkSize, numBS));
    bs_phasorSumIso = complex(zeros(chunkSize, numBS));

    for j = 1:chunkSize
        raysForRx = raysAllTxRx(:, j);
        hasRays = ~cellfun(@isempty, raysForRx);

        if ~any(hasRays)
            continue;
        end

        raysValid = raysForRx(hasRays);
        txValid   = txSite(hasRays);
        gainValid = gainTx(hasRays);

        if isscalar(allocatedTxPower)
            allocPowerValid = allocatedTxPower;
        else
            allocPowerValid = allocatedTxPower(hasRays);
        end

        exp = calcPower4(txValid, rxChunk(j), raysValid, fc, allocPowerValid, gainValid, gainRx);

        % Totals
        tot_incPow(j) = exp.totalIncoherentPower_dBm;
        tot_cohPow(j) = exp.totalCoherentPower_dBm;
        tot_incEf(j)  = exp.totalIncoherentEfield_dBuv;
        tot_cohEf(j)  = exp.totalCoherentEfield_dBuv;

        % Per-BS — map back to original BS indices
        bsIdx = find(hasRays);
        for b = 1:numel(bsIdx)
            bs_cohPow(j, bsIdx(b))    = exp.coherentPower_dBm(b);
            bs_incPow(j, bsIdx(b))    = exp.incoherentPower_dBm(b);
            bs_cohPowIso(j, bsIdx(b)) = exp.coherentPowerIso_dBm(b);
            bs_incPowIso(j, bsIdx(b)) = exp.incoherentPowerIso_dBm(b);
            bs_cohEf(j, bsIdx(b))     = exp.coherentEfield_dBuv(b);
            bs_incEf(j, bsIdx(b))     = exp.incoherentEfield_dBuv(b);
            bs_phasorSum(j, bsIdx(b))    = exp.phasorSum_perBS(b);
            bs_phasorSumIso(j, bsIdx(b)) = exp.phasorSumIso_perBS(b);
        end
    end

    % Save block
    m.totalIncoherentPower_dBm(idxStart:idxEnd, 1) = tot_incPow;
    m.totalCoherentPower_dBm(idxStart:idxEnd, 1)   = tot_cohPow;
    m.totalIncoherentEfield_dBuv(idxStart:idxEnd, 1) = tot_incEf;
    m.totalCoherentEfield_dBuv(idxStart:idxEnd, 1)   = tot_cohEf;

    m.coherentPower_dBm(idxStart:idxEnd, 1:numBS)      = bs_cohPow;
    m.incoherentPower_dBm(idxStart:idxEnd, 1:numBS)    = bs_incPow;
    m.coherentPowerIso_dBm(idxStart:idxEnd, 1:numBS)   = bs_cohPowIso;
    m.incoherentPowerIso_dBm(idxStart:idxEnd, 1:numBS) = bs_incPowIso;
    m.coherentEfield_dBuv(idxStart:idxEnd, 1:numBS)    = bs_cohEf;
    m.incoherentEfield_dBuv(idxStart:idxEnd, 1:numBS)  = bs_incEf;

    m.phasorSum_perBS(idxStart:idxEnd, 1:numBS)    = bs_phasorSum;
    m.phasorSumIso_perBS(idxStart:idxEnd, 1:numBS) = bs_phasorSumIso;

    % Progress
    blockTime = toc(tBlock);
    elapsed = toc(tStart);
    pctDone = idxEnd / N * 100;
    remaining = (elapsed / idxEnd) * (N - idxEnd);
    fprintf('Block %d/%d | %d-%d of %d (%.1f%%) | Block: %.1fs | Elapsed: %.1fm | ETA: %.1fm\n', ...
        blockNum, totalBlocks, idxStart, idxEnd, N, pctDone, ...
        blockTime, elapsed/60, remaining/60);
end

% Build output —> totals
powerMatrix.totalIncoherentPower_dBm = m.totalIncoherentPower_dBm;
powerMatrix.totalCoherentPower_dBm   = m.totalCoherentPower_dBm;
powerMatrix.totalIncoherentEfield_dBuv = m.totalIncoherentEfield_dBuv;
powerMatrix.totalCoherentEfield_dBuv   = m.totalCoherentEfield_dBuv;

% Build output —> per BS
powerMatrix.coherentPower_dBm      = m.coherentPower_dBm;
powerMatrix.incoherentPower_dBm    = m.incoherentPower_dBm;
powerMatrix.coherentPowerIso_dBm   = m.coherentPowerIso_dBm;
powerMatrix.incoherentPowerIso_dBm = m.incoherentPowerIso_dBm;
powerMatrix.coherentEfield_dBuv    = m.coherentEfield_dBuv;
powerMatrix.incoherentEfield_dBuv  = m.incoherentEfield_dBuv;

powerMatrix.phasorSum_perBS    = m.phasorSum_perBS;
powerMatrix.phasorSumIso_perBS = m.phasorSumIso_perBS;

% Coordinates
powerMatrix.lat = lat;
powerMatrix.lon = lon;

fprintf('Power calculation complete.\n');
end
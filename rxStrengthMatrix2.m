% function [powerMatrix] = rxStrengthMatrix2(rxSamplingArray, txSite, pm, fc, allocatedTxPower)
% % RXSTRENGTHMATRIX Calculates received power/E-field at each sampling point
% %   Wraps calcPower in a chunked loop for memory efficiency.
% 
% N = numel(rxSamplingArray);
% lat = [rxSamplingArray.Latitude]';
% lon = [rxSamplingArray.Longitude]';
% 
% % Setup incremental save file
% resultFolder = 'C:\Users\chris\OneDrive - Students RWTH Aachen University\RWTH\Master Thesis\Results';
% saveFile = fullfile(resultFolder, 'power_progress.mat');
% m = matfile(saveFile, 'Writable', true);
% m.incoherentPower_dBm = -inf(N, 1);
% m.coherentPower_dBm   = -inf(N, 1);
% m.incoherentEfield_dBuv = -inf(N, 1);
% m.coherentEfield_dBuv   = -inf(N, 1);
% 
% blockSize = 100;
% fprintf('Starting power calculation for %d sites...\n', N);
% 
% for i = 1:blockSize:N
% 
%     idxStart = i;
%     idxEnd = min(i + blockSize - 1, N);
%     chunkSize = idxEnd - idxStart + 1;
%     rxChunk = rxSamplingArray(idxStart:idxEnd).';
% 
%     % Ray trace this chunk: single Tx → multiple Rx
%     raysPerChunk = raytrace(txSite, rxChunk, pm);
% 
%     % Pre-allocate chunk results
%     chunkIncPow   = -inf(chunkSize, 1);
%     chunkCohPow   = -inf(chunkSize, 1);
%     chunkIncEf    = -inf(chunkSize, 1);
%     chunkCohEf    = -inf(chunkSize, 1);
% 
%     for j = 1:chunkSize
% 
%         raysForThisRx = raysPerChunk(:, j);
% 
%         if isempty(raysPerChunk{j})
%             continue;  % stays -Inf
%         end
% 
%         % Wrap single-Tx rays into cell for calcPower
%         % calcPower expects rays{i} for each Tx, so pass {rays} with one element
%         singleRxExposure = calcPower(txSite, rxChunk(j), raysForThisRx, fc, allocatedTxPower);
% 
%         chunkIncPow(j) = singleRxExposure.totalIncoherentPower_dBm;
%         chunkCohPow(j) = singleRxExposure.totalCoherentPower_dBm;
%         chunkIncEf(j)  = singleRxExposure.totalIncoherentEfield_dBuv;
%         chunkCohEf(j)  = singleRxExposure.totalCoherentEfield_dBuv;
%     end
% 
%     % Save to disk incrementally
%     m.incoherentPower_dBm(idxStart:idxEnd, 1) = chunkIncPow;
%     m.coherentPower_dBm(idxStart:idxEnd, 1) = chunkCohPow;
%     m.incoherentEfield_dBuv(idxStart:idxEnd, 1) = chunkIncEf;
%     m.coherentEfield_dBuv(idxStart:idxEnd, 1)= chunkCohEf;
% 
%     clear rxChunk raysPerChunk;
%     fprintf('Processed %d to %d (%.1f%%)\n', idxStart, idxEnd, (idxEnd/N)*100);
% end
% 
% % Build output struct
% powerMatrix.incoherentPower_dBm   = m.incoherentPower_dBm;
% powerMatrix.coherentPower_dBm     = m.coherentPower_dBm;
% powerMatrix.incoherentEfield_dBuv = m.incoherentEfield_dBuv;
% powerMatrix.coherentEfield_dBuv   = m.coherentEfield_dBuv;
% powerMatrix.lat = lat;
% powerMatrix.lon = lon;
% 
% fprintf('Power calculation complete.\n');
% end

% function [powerMatrix] = rxStrengthMatrix2(rxSamplingArray, txSite, pm, fc, allocatedTxPower)
% 
% N = numel(rxSamplingArray);
% lat = [rxSamplingArray.Latitude]';
% lon = [rxSamplingArray.Longitude]';
% 
% % Pre-allocate results
% incPow  = -inf(N, 1);
% cohPow  = -inf(N, 1);
% incEf   = -inf(N, 1);
% cohEf   = -inf(N, 1);
% 
% fprintf('Starting power calculation for %d points...\n', N);
% 
% parfor i = 1:N
%     rx_i = rxSamplingArray(i);
% 
%     % Ray trace: all Tx to this one Rx
%     raysForRx = raytrace(txSite, rx_i, pm);  % {numTx × 1} cell
% 
%     hasRays = ~cellfun(@isempty, raysForRx);
%     if ~any(hasRays)
%         continue;
%     end
% 
%     raysValid = raysForRx(hasRays);
%     txValid   = txSite(hasRays);
% 
%     exp = calcPower(txValid, rx_i, raysValid, fc, allocatedTxPower);
% 
%     incPow(i) = exp.totalIncoherentPower_dBm;
%     cohPow(i) = exp.totalCoherentPower_dBm;
%     incEf(i)  = exp.totalIncoherentEfield_dBuv;
%     cohEf(i)  = exp.totalCoherentEfield_dBuv;
% end
% 
% % Build output
% powerMatrix.incoherentPower_dBm   = incPow;
% powerMatrix.coherentPower_dBm     = cohPow;
% powerMatrix.incoherentEfield_dBuv = incEf;
% powerMatrix.coherentEfield_dBuv   = cohEf;
% powerMatrix.lat = lat;
% powerMatrix.lon = lon;
% 
% fprintf('Power calculation complete.\n');
% end

function [powerMatrix] = rxStrengthMatrix2(rxSamplingArray, txSite, pm, fc, allocatedTxPower)

N = numel(rxSamplingArray);
lat = [rxSamplingArray.Latitude]';
lon = [rxSamplingArray.Longitude]';

% Setup incremental save
resultFolder = 'C:\Users\chris\OneDrive - Students RWTH Aachen University\RWTH\Master Thesis\Results';
saveFile = fullfile(resultFolder, 'power_progress.mat');
m = matfile(saveFile, 'Writable', true);
m.incoherentPower_dBm   = -inf(N, 1);
m.coherentPower_dBm     = -inf(N, 1);
m.incoherentEfield_dBuv = -inf(N, 1);
m.coherentEfield_dBuv   = -inf(N, 1);

blockSize = 200;  % tune this: big enough for parfor efficiency, small enough for frequent saves
fprintf('Starting power calculation for %d points...\n', N);

for block = 1:blockSize:N
    idxStart = block;
    idxEnd   = min(block + blockSize - 1, N);
    chunkSize = idxEnd - idxStart + 1;

    rxChunk = rxSamplingArray(idxStart:idxEnd);

    % Pre-allocate chunk results
    incPow = -inf(chunkSize, 1);
    cohPow = -inf(chunkSize, 1);
    incEf  = -inf(chunkSize, 1);
    cohEf  = -inf(chunkSize, 1);

    parfor j = 1:chunkSize
        rx_j = rxChunk(j);
        raysForRx = raytrace(txSite, rx_j, pm);  % {numTx × 1}

        hasRays = ~cellfun(@isempty, raysForRx);
        if ~any(hasRays)
            continue;
        end

        raysValid = raysForRx(hasRays);
        txValid   = txSite(hasRays);

        exp = calcPower(txValid, rx_j, raysValid, fc, allocatedTxPower);

        incPow(j) = exp.totalIncoherentPower_dBm;
        cohPow(j) = exp.totalCoherentPower_dBm;
        incEf(j)  = exp.totalIncoherentEfield_dBuv;
        cohEf(j)  = exp.totalCoherentEfield_dBuv;
    end

    % Save block to disk (outside parfor)
    m.incoherentPower_dBm(idxStart:idxEnd, 1)   = incPow;
    m.coherentPower_dBm(idxStart:idxEnd, 1)     = cohPow;
    m.incoherentEfield_dBuv(idxStart:idxEnd, 1) = incEf;
    m.coherentEfield_dBuv(idxStart:idxEnd, 1)   = cohEf;

    fprintf('Processed %d to %d of %d (%.1f%%)\n', idxStart, idxEnd, N, (idxEnd/N)*100);
end

% Build output
powerMatrix.incoherentPower_dBm   = m.incoherentPower_dBm;
powerMatrix.coherentPower_dBm     = m.coherentPower_dBm;
powerMatrix.incoherentEfield_dBuv = m.incoherentEfield_dBuv;
powerMatrix.coherentEfield_dBuv   = m.coherentEfield_dBuv;
powerMatrix.lat = lat;
powerMatrix.lon = lon;

fprintf('Power calculation complete.\n');
end
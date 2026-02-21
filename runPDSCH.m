function [throughputMbps] = runPDSCH(channelObj, numFrames, numLayer, numHARQ, SNRdB, usedMod, numRB, coderate, SCS, weightsTx, numRBalloc)
%%RUNPDSCH calculates download speed
% based on https://de.mathworks.com/help/5g/ug/nr-pdsch-throughput.html
% Most of the settings are standard and taken directly from the example in
% the link
tic
release(channelObj) % Enables configuration of nrCDLchannel system object

channelObj.ChannelResponseOutput = 'ofdm-response'; % Enables Rx channel estimation aligned to OFDM symbols

% Dont give additional diagnostics?
DisplayDiagnostics = true; % true/false
DisplaySimulationInformation = true;

% Try GPU acceleration
% UseGPU = "on";

%% Define OFDM parameter and resources / Carrier and PDSCH Configuration
% Instead of prod(channelObj.TransmitAntennaArray.Size) takes polarization
% into account

% Calculate the number of spatial positions (geometry)
numSpatial = prod(channelObj.TransmitAntennaArray.Size);

% Calculate the number of polarizations (based on your ElementSet)
numPol = numel(channelObj.TransmitAntennaArray.ElementSet);

% Calculate total Transmit Antenna Ports (Logic Ports)
NumTxAnts = numSpatial * numPol; % This will yield 128

% NumTxAnts = prod(channelObj.TransmitAntennaArray.Size);
% NumRxAnts = prod(channelObj.ReceiveAntennaArray.Size);

% Set waveform type
carrier = nrCarrierConfig; % Creates 5G NR carrier configuration object
carrier.NSizeGrid = numRB; % BW in number of resource blocks
carrier.SubcarrierSpacing = SCS; % 30 kHz
carrier.CyclicPrefix = 'Normal'; % if need change CP, 'Normal' or 'Extended' (Extended CP is relevant for 60 kHz SCS only)
carrier.NCellID = 1; % Cell identity

% Creates PDSCH configuration and an additional struct for extra parameters
% This PDSCH definition is the basis for all PDSCH transmissions in the BLER (Block Error Rate) simulation
pdsch = nrPDSCHConfig;
pdschExtension = struct(); % This structure is to hold additional simulation parameters for the DL-SCH and PDSC

%% Set allocated RBs
% Define PDSCH time-frequency resource allocation per slot to be full grid
pdsch.PRBSet = 0:(numRBalloc-1); % Full frequencies range, PDSCH PRB allocation
pdsch.SymbolAllocation = [0 carrier.SymbolsPerSlot]; % PDSCH uses all OFDM symbols per slot
pdsch.MappingType = 'A'; % PDSCH mapping type ('A'(slot-wise),'B'(non slot-wise))


% PDSCH resource block mapping (TS 38.211 Section 7.3.1.6)
pdsch.VRBToPRBInterleaving = 0; % Disable interleaved resource mapping
pdsch.VRBBundleSize = 4;

% Define the number of transmission layers to be used
pdsch.NumLayers = numLayer; % Number of PDSCH transmission layers

% Set used modulation and coderate based on SNR and MCS index
% pdsch.Modulation = usedMod;
% pdschExtension.TargetCodeRate = coderate;

% Define codeword modulation and target coding rate
% The number of codewords is directly dependent on the number of layers so ensure that
% layers are set first before getting the codeword number

if pdsch.NumCodewords > 1                             % Multicodeword transmission (when number of layers being > 4)
    pdsch.Modulation = {usedMod,usedMod};             % 'QPSK', '16QAM', '64QAM', '256QAM'
    pdschExtension.TargetCodeRate = [coderate coderate];   % Code rate used to calculate transport block sizes
else
    pdsch.Modulation = usedMod;                       % 'QPSK', '16QAM', '64QAM', '256QAM'
    pdschExtension.TargetCodeRate = coderate;         % Code rate used to calculate transport block sizes
end

% DM-RS and antenna port configuration (TS 38.211 Section 7.4.1.1)
% DeModulation Reference Signal
pdsch.DMRS.DMRSPortSet = 0:pdsch.NumLayers-1; % DM-RS ports to use for the layers
pdsch.DMRS.DMRSTypeAPosition = 2; % Mapping type A only. First DM-RS symbol position (2,3)
pdsch.DMRS.DMRSConfigurationType = 1; % DM-RS configuration type (1,2)
pdsch.DMRS.DMRSLength = 1; % Number of front-loaded DM-RS symbols (1(single symbol),2(double symbol)) -> use 2 to be able to use 8 layers
pdsch.DMRS.DMRSAdditionalPosition = 2; % Additional DM-RS symbol positions (max range 0...3) length 2 and positions 2 doesnt work
pdsch.DMRS.NumCDMGroupsWithoutData = 2; % Number of CDM groups without data

% PT-RS configuration (TS 38.211 Section 7.4.1.2)
% Phase Tracking Reference Signal
pdsch.EnablePTRS = 0; % Enable or disable PT-RS (1 or 0)
pdsch.PTRS.TimeDensity = 1; % PT-RS time density (L_PT-RS) (1, 2, 4)
pdsch.PTRS.FrequencyDensity = 2; % PT-RS frequency density (K_PT-RS) (2 or 4)
pdsch.PTRS.REOffset = '00'; % PT-RS resource element offset ('00', '01', '10', '11')
pdsch.PTRS.PTRSPortSet = []; % PT-RS antenna port, subset of DM-RS port set. Empty corresponds to lower DM-RS port number

% Reserved PRB patterns, if required (for CORESETs, forward compatibility etc)
pdsch.ReservedPRB{1}.SymbolSet = []; % Reserved PDSCH symbols
pdsch.ReservedPRB{1}.PRBSet = []; % Reserved PDSCH PRBs
pdsch.ReservedPRB{1}.Period = []; % Periodicity of reserved resources

%% Additional simulation and DL-SCH related parameters
% PDSCH PRB bundling (TS 38.214 Section 5.1.2.3)
pdschExtension.PRGBundleSize = []; % 2, 4, or [] to signify "wideband"

% HARQ process and rate matching/TBS parameters
pdschExtension.XOverhead = 6 * pdsch.EnablePTRS; % Set PDSCH rate matching overhead for TBS (Xoh) to 6 when PT-RS is enabled, otherwise 0
pdschExtension.NHARQProcesses = numHARQ; % Number of parallel HARQ processes to use
pdschExtension.EnableHARQ = true; % Enable retransmissions for each process, using RV sequence [0,2,3,1]

% LDPC decoder parameters
% Available algorithms: 'Belief propagation', 'Layered belief propagation', 'Normalized min-sum', 'Offset min-sum'
pdschExtension.LDPCDecodingAlgorithm = 'Normalized min-sum';
pdschExtension.MaximumLDPCIterationCount = 6;

% Define data type ('single' or 'double') for resource grids and waveforms
DataType = 'single';

% Get information about the baseband waveform after OFDM modulation step
waveformInfo = nrOFDMInfo(carrier);

% Get channel information and max channel delay
chInfo = info(channelObj);
maxChDelay = chInfo.MaximumChannelDelay;

%% Processing Loop

% Array to store the maximum throughput for all SNR points
maxThroughput = zeros(1,1);
% Array to store the simulation throughput for all SNR points
simThroughput = zeros(1,1);

% Set up redundancy version (RV) sequence for all HARQ processes
if pdschExtension.EnableHARQ
    % In the final report of RAN WG1 meeting #91 (R1-1719301), it was
    % observed in R1-1717405 that if performance is the priority, [0 2 3 1]
    % should be used. If self-decodability is the priority, it should be
    % taken into account that the upper limit of the code rate at which
    % each RV is self-decodable is in the following order: 0>3>2>1
    rvSeq = [0 2 3 1];
else
    % HARQ disabled - single transmission with RV=0, no retransmissions
    rvSeq = 0;
end

% Create DL-SCH encoder system object to perform transport channel encoding
encodeDLSCH = nrDLSCH;
encodeDLSCH.MultipleHARQProcesses = true;

% MATLAB wants coderate < 1, table gives them as integers -> coderate/1024
% provided in function handover
pdschExtension.TargetCodeRate = pdschExtension.TargetCodeRate;
encodeDLSCH.TargetCodeRate = pdschExtension.TargetCodeRate;

% Create DL-SCH decoder system object to perform transport channel decoding
% Use layered belief propagation for LDPC decoding, with half the number of
% iterations as compared to the default for belief propagation decoding
decodeDLSCH = nrDLSCHDecoder;
decodeDLSCH.MultipleHARQProcesses = true;
decodeDLSCH.TargetCodeRate = pdschExtension.TargetCodeRate;
decodeDLSCH.LDPCDecodingAlgorithm = pdschExtension.LDPCDecodingAlgorithm;
decodeDLSCH.MaximumLDPCIterationCount = pdschExtension.MaximumLDPCIterationCount;

% Reset the random number generator so that each SNR point will
% experience the same noise realization
rng('default');

pdschextra = pdschExtension;
pathFilters = [];

% Specify fixed order in which to cycle through the HARQ process IDs
harqSequence = 0:pdschextra.NHARQProcesses-1;

% Initialize the state of all HARQ processes
harqEntity = HARQEntity(harqSequence ,rvSeq, pdsch.NumCodewords);

% Reset the channel so that each SNR point will experience the same
% channel realization
reset(channelObj);

% Total number of slots in the simulation period
NSlots = numFrames * carrier.SlotsPerFrame;

% Timing offset, updated in every slot for perfect synchronization and
% when the correlation is strong for practical synchronization
offset = 0;


% % 1. Measure the Transmit Power (Average power per antenna)
% % We use the Tx waveform because it is stable and clean.
% txSigPower = mean(abs(txWaveform(:)).^2);
%
%
% % Noise power, normalized by the IFFT size used in OFDM modulation, as
% % the OFDM modulator applies this normalization to the transmitted
% % waveform. Also normalize by the number of receive antennas, as the
% % channel model applies this normalization to the received waveform by
% % default. Calculate the noise power per RE to act as the noise
% % estimate if perfect channel estimation is enabled
SNR_linear = 10^(SNRdB/10);
% noisePower = txSigPower / SNR_linear;
% noiseStd = sqrt(noisePower/2);

% N0 = 1/sqrt(NumRxAnts * double(waveformInfo.Nfft) * SNR_linear);
% % nVar = N0^2*double(waveformInfo.Nfft);

% Loop over the entire waveform length
% NSlots depends on numFrames

%estChannelGridAnts = getInitialChannelEstimate(carrier, channelObj, DataType, maxChDelay);
%weightsTx = hSVDPrecoders(carrier, pdsch, estChannelGridAnts, pdschExtension.PRGBundleSize);


for nslot = 0:NSlots-1

    fprintf("Progress: Slot %d from %d \n", nslot, NSlots)
    % Update the carrier slot numbers for new slot
    carrier.NSlot = nslot;

    % Calculate the transport block sizes for the transmission in the slot
    [pdschIndices, pdschIndicesInfo] = nrPDSCHIndices(carrier, pdsch);
    trBlkSizes = nrTBS(usedMod, numLayer, numel(pdsch.PRBSet),pdschIndicesInfo.NREPerPRB, pdschExtension.TargetCodeRate, pdschExtension.XOverhead);

    % HARQ processing
    for cwIdx = 1:pdsch.NumCodewords

        % If new data for current process and codeword then create a new DL-SCH transport block
        if harqEntity.NewData(cwIdx)
            trBlk = randi([0 1], trBlkSizes(cwIdx), 1);
            setTransportBlock(encodeDLSCH, trBlk, cwIdx - 1, harqEntity.HARQProcessID);

            % If new data because of previous RV sequence time out then flush decoder soft buffer explicitly
            if harqEntity.SequenceTimeout(cwIdx)
                resetSoftBuffer(decodeDLSCH, cwIdx - 1, harqEntity.HARQProcessID);
            end
        end
    end

    % Encode the DL-SCH transport blocks
    codedTrBlocks = encodeDLSCH(pdsch.Modulation, pdsch.NumLayers,...
        pdschIndicesInfo.G, harqEntity.RedundancyVersion, harqEntity.HARQProcessID);

    % Create a resource grid like the prototype array for a slot
    pdschGrid = nrResourceGrid(carrier, NumTxAnts, OutputDataType = DataType);

    % PDSCH modulation and precoding
    pdschSymbols = nrPDSCH(carrier, pdsch, codedTrBlocks);
    [pdschAntSymbols, pdschAntIndices] = nrPDSCHPrecode(carrier, pdschSymbols, pdschIndices, weightsTx);

    % PDSCH mapping in grid associated with PDSCH transmission period
    pdschGrid(pdschAntIndices) = pdschAntSymbols;

    % PDSCH DM-RS precoding and mapping
    dmrsSymbols = nrPDSCHDMRS(carrier, pdsch);
    dmrsIndices = nrPDSCHDMRSIndices(carrier, pdsch);
    [dmrsAntSymbols, dmrsAntIndices] = nrPDSCHPrecode(carrier, dmrsSymbols, dmrsIndices, weightsTx);
    pdschGrid(dmrsAntIndices) = dmrsAntSymbols;

    % PDSCH PT-RS precoding and mapping
    ptrsSymbols = nrPDSCHPTRS(carrier, pdsch);
    ptrsIndices = nrPDSCHPTRSIndices(carrier, pdsch);
    [ptrsAntSymbols, ptrsAntIndices] = nrPDSCHPrecode(carrier, ptrsSymbols, ptrsIndices, weightsTx);
    pdschGrid(ptrsAntIndices) = ptrsAntSymbols;

    % OFDM modulation
    txWaveform = nrOFDMModulate(carrier, pdschGrid);

    % Pass data through channel model. Append zeros at the end of the
    % transmitted waveform to flush channel content. These zeros take
    % into account any delay introduced in the channel. This is a mix
    % of multipath delay and implementation delay. This value may
    % change depending on the sampling rate, delay profile, and delay
    % spread. The channel model also returns the OFDM channel response
    % and timing offset for the specified carrier
    txWaveform = [txWaveform; zeros(maxChDelay, size(txWaveform, 2))];
    [rxWaveform, ofdmResponse, timingOffset] = channelObj(txWaveform, carrier);
    % [rxWaveform, ofdmResponse, timingOffset] = channelObj(txWaveform);
    % Add AWGN to the received time domain waveform

    % sigPower = mean(abs(rxWaveform(:)).^2);
    % %%% 2) Noise-Power = signalPower / SNR
    % noisePower = sigPower / SNR_linear;
    % %%% 3) sqrt(NoisePower/2) => da komplexes AWGN (I + Q)
    % noiseStd = sqrt(noisePower/2);
    % %%% 4) Erzeuge das AWGN
    % noise = noiseStd * (randn(size(rxWaveform)) + 1i*randn(size(rxWaveform)));
    %
    %
    %
    % %%% 5) Addiere auf das empfangene Signal
    % rxWaveform = rxWaveform + noise;
    % Why not these? Not imaginary?
    % noise = N0 * (randn(size(rxWaveform)) + 1i*randn(size(rxWaveform)));
    %
    % rxWaveform = rxWaveform + noise;

    rxSig = rxWaveform(1:end-maxChDelay, :);
    sigPower = mean(abs(rxSig(:)).^2);

    % Noise power relative to actual received signal
    noisePower = sigPower / SNR_linear;

    % For REAL-valued noise (matching the reference example convention):
    noise = sqrt(noisePower) * randn(size(rxWaveform), 'like', rxWaveform);
    rxWaveform = rxWaveform + noise;


    % Practical synchronization. Correlate the received waveform
    % with the PDSCH DM-RS to give timing offset estimate 't' and
    % correlation magnitude 'mag'. The function
    % hSkipWeakTimingOffset is used to update the receiver timing
    % offset. If the correlation peak in 'mag' is weak, the current
    % timing estimate 't' is ignored and the previous estimate
    % 'offset' is used
    [t, mag] = nrTimingEstimate(carrier, rxWaveform, dmrsIndices, dmrsSymbols);
    offset = hSkipWeakTimingOffset(offset, t, mag);

    % Display a warning if the estimated timing offset exceeds the
    % maximum channel delay
    if offset > maxChDelay
        warning(['Estimated timing offset (%d) is greater than the maximum channel delay (%d).' ...
            ' This will result in a decoding failure. This may be caused by low SNR,' ...
            ' or not enough DM-RS symbols to synchronize successfully.'],offset,maxChDelay);
    end

    rxWaveform = rxWaveform(1+offset:end,:);

    % Perform OFDM demodulation on the received data to recreate the
    % resource grid, including padding in the event that practical
    % synchronization results in an incomplete slot being demodulated
    rxGrid = nrOFDMDemodulate(carrier, rxWaveform);
    [K, L, R] = size(rxGrid);
    if (L < carrier.SymbolsPerSlot)
        rxGrid = cat(2, rxGrid, zeros(K, carrier.SymbolsPerSlot - L, R));
    end

    % Practical channel estimation between the received grid and
    % each transmission layer, using the PDSCH DM-RS for each
    % layer. This channel estimate includes the effect of
    % transmitter precoding
    [estChannelGridPorts, noiseEst] = hSubbandChannelEstimate(carrier, rxGrid, dmrsIndices, dmrsSymbols, pdschextra.PRGBundleSize, 'CDMLengths', pdsch.DMRS.CDMLengths);

    % Average noise estimate across PRGs and layers
    noiseEst = mean(noiseEst,'all');

    % Get PDSCH resource elements from the received grid and
    % channel estimate
    [pdschRx, pdschHest] = nrExtractResources(pdschIndices, rxGrid, estChannelGridPorts);

    % Remove precoding from estChannelGridPorts to get channel
    % estimate w.r.t. antennas
    estChannelGridAnts = precodeChannelEstimate(carrier, estChannelGridPorts, conj(weightsTx)); % Why transpose?

    % Diagnostic: check channel spatial condition
    H_mid = squeeze(estChannelGridAnts(round(size(estChannelGridAnts,1)/2), 1, :, :));
    fprintf('H_mid size: %s\n', mat2str(size(H_mid)));
    if ismatrix(H_mid)
        sv = svd(H_mid);
        fprintf('Singular values: ');
        fprintf('%.6f  ', sv);
        fprintf('\nCondition (sv1/sv_end): %.1f\n', sv(1)/sv(end));
    end




    % Equalization
    % Look how the csi is derived!
    [pdschEq, csi] = nrEqualizeMMSE(pdschRx, pdschHest, noiseEst);

    % Common phase error (CPE) compensation
    if ~isempty(ptrsIndices)

        % Initialize temporary grid to store equalized symbols
        tempGrid = nrResourceGrid(carrier, pdsch.NumLayers);

        % Extract PT-RS symbols from received grid and estimated
        % channel grid
        [ptrsRx, ptrsHest, ~, ~, ptrsHestIndices, ptrsLayerIndices] = nrExtractResources(ptrsIndices, rxGrid, estChannelGridAnts, tempGrid);
        ptrsHest = nrPDSCHPrecode(carrier, ptrsHest, ptrsHestIndices, permute(weightsTx, [2 1 3]));

        % Equalize PT-RS symbols and map them to tempGrid
        ptrsEq = nrEqualizeMMSE(ptrsRx, ptrsHest, noiseEst);
        tempGrid(ptrsLayerIndices) = ptrsEq;

        % Estimate the residual channel at the PT-RS locations in
        % tempGrid
        cpe = nrChannelEstimate(tempGrid, ptrsIndices, ptrsSymbols);

        % Sum estimates across subcarriers, receive antennas, and
        % layers. Then, get the CPE by taking the angle of the
        % resultant sum
        cpe = angle(sum(cpe, [1 3 4]));

        % Map the equalized PDSCH symbols to tempGrid
        tempGrid(pdschIndices) = pdschEq;

        % Correct CPE in each OFDM symbol within the range of reference
        % PT-RS OFDM symbols
        symLoc = pdschIndicesInfo.PTRSSymbolSet(1) + 1:pdschIndicesInfo.PTRSSymbolSet(end) + 1;
        tempGrid(:, symLoc, :) = tempGrid(:, symLoc, :).*exp(-1i*cpe(symLoc));

        % Extract PDSCH symbols
        pdschEq = tempGrid(pdschIndices);
    end

    % Decode PDSCH physical channel
    [dlschLLRs, rxSymbols] = nrPDSCHDecode(carrier, pdsch, pdschEq, noiseEst);

    % Display EVM per layer, per slot and per RB
    if (DisplayDiagnostics)
        plotLayerEVM(NSlots, nslot, pdsch, size(pdschGrid), pdschIndices, pdschSymbols, pdschEq);
    end

    % Scale LLRs by CSI
    csi = nrLayerDemap(csi); % CSI layer demapping
    for cwIdx = 1:pdsch.NumCodewords
        Qm = length(dlschLLRs{cwIdx}) / length(rxSymbols{cwIdx}); % bits per symbol
        csi{cwIdx} = repmat(csi{cwIdx}.', Qm, 1); % expand by each bit per symbol
        dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:); % scale by CSI
    end

    % Decode the DL-SCH transport channel
    decodeDLSCH.TransportBlockLength = trBlkSizes;
    [decbits, blkerr] = decodeDLSCH(dlschLLRs, pdsch.Modulation, pdsch.NumLayers, harqEntity.RedundancyVersion, harqEntity.HARQProcessID);

    % Store values to calculate throughput
    simThroughput(1) = simThroughput(1) + sum(~blkerr .* trBlkSizes);
    maxThroughput(1) = maxThroughput(1) + sum(trBlkSizes);

    % Update current process with CRC error and advance to next process
    procstatus = updateAndAdvance(harqEntity, blkerr, trBlkSizes, pdschIndicesInfo.G);
    if (DisplaySimulationInformation)
        fprintf('\n(%3.2f%%) NSlot=%d, %s', 100*(nslot + 1) / NSlots, nslot, procstatus);
    end

    % % Get precoding matrix for next slot
    newWtx = hSVDPrecoders(carrier, pdsch, estChannelGridAnts, pdschextra.PRGBundleSize);
    weightsTx = newWtx;



end

% Display the results dynamically in the command window
if (DisplaySimulationInformation)
    fprintf('\n');
end
fprintf('Throughput (Mbps) for %d frame(s) = %.4f\n', numFrames, 1e-6 * simThroughput / (numFrames * 10e-3));
fprintf('Throughput (%%) for %d frame(s) = %.4f\n\n', numFrames, simThroughput * 100 / maxThroughput);
throughputMbps = 1e-6 * simThroughput / (numFrames * 10e-3);

toc
end



%% Local functions
function validateNumLayers(simParameters)
% Validate the number of layers, relative to the antenna geometry

numlayers = simParameters.PDSCH.NumLayers;
ntxants = simParameters.NTxAnts;
nrxants = simParameters.NRxAnts;
antennaDescription = sprintf('min(NTxAnts,NRxAnts) = min(%d,%d) = %d',ntxants,nrxants,min(ntxants,nrxants));
if numlayers > min(ntxants,nrxants)
    error('The number of layers (%d) must satisfy NumLayers <= %s', ...
        numlayers,antennaDescription);
end

% Display a warning if the maximum possible rank of the channel equals
% the number of layers
if (numlayers > 2) && (numlayers == min(ntxants,nrxants))
    warning(['The maximum possible rank of the channel, given by %s, is equal to NumLayers (%d).' ...
        ' This may result in a decoding failure under some channel conditions.' ...
        ' Try decreasing the number of layers or increasing the channel rank' ...
        ' (use more transmit or receive antennas).'],antennaDescription,numlayers); %#ok<SPWRN>
end

end

function estChannelGrid = getInitialChannelEstimate(carrier,propchannel,dataType,maxChDelay)
% Obtain channel estimate before first transmission. This can be used to
% obtain a precoding matrix for the first slot.

ofdmInfo = nrOFDMInfo(carrier);

% Clone of the channel
chClone = propchannel.clone();
chClone.release();

% No filtering needed to get perfect channel estimate
chClone.ChannelFiltering = false;
chClone.OutputDataType = dataType;
chClone.NumTimeSamples = (ofdmInfo.SampleRate/1000/carrier.SlotsPerSubframe)+maxChDelay;

% Get the perfect channel estimate
estChannelGrid = chClone(carrier);

end

function estChannelGrid = precodeChannelEstimate(carrier,estChannelGrid,W)
% Apply precoding matrix W to the last dimension of the channel estimate

[K,L,R,P] = size(estChannelGrid);
estChannelGrid = reshape(estChannelGrid,[K*L R P]);
estChannelGrid = nrPDSCHPrecode(carrier,estChannelGrid,reshape(1:numel(estChannelGrid),[K*L R P]),W);
estChannelGrid = reshape(estChannelGrid,K,L,R,[]);

end

function plotLayerEVM(NSlots,nslot,pdsch,siz,pdschIndices,pdschSymbols,pdschEq)
% Plot EVM information

persistent slotEVM;
persistent rbEVM
persistent evmPerSlot;

if (nslot==0)
    slotEVM = comm.EVM;
    rbEVM = comm.EVM;
    evmPerSlot = NaN(NSlots,pdsch.NumLayers);
    figure;
end
evmPerSlot(nslot+1,:) = slotEVM(pdschSymbols,pdschEq);
subplot(2,1,1);
plot(0:(NSlots-1),evmPerSlot,'o-');
xlabel('Slot number');
ylabel('EVM (%)');
legend("layer " + (1:pdsch.NumLayers),'Location','EastOutside');
title('EVM per layer per slot');

subplot(2,1,2);
[k,~,p] = ind2sub(siz,pdschIndices);
rbsubs = floor((k-1) / 12);
NRB = siz(1) / 12;
evmPerRB = NaN(NRB,pdsch.NumLayers);
for nu = 1:pdsch.NumLayers
    for rb = unique(rbsubs).'
        this = (rbsubs==rb & p==nu);
        evmPerRB(rb+1,nu) = rbEVM(pdschSymbols(this),pdschEq(this));
    end
end
plot(0:(NRB-1),evmPerRB,'x-');
xlabel('Resource block');
ylabel('EVM (%)');
legend("layer " + (1:pdsch.NumLayers),'Location','EastOutside');
title(['EVM per layer per resource block, slot #' num2str(nslot)]);

drawnow;

end



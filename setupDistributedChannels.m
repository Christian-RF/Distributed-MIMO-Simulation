function [channelData, hEstCollective] = setupDistributedChannels(rays, tx, rx, fc, UEOrientation, numRB, SCS, numSlot)
%SETUPDISTRIBUTEDCHANNELS Sets up CDL channels and calculates perfect estimates
%   Returns a struct containing channel objects and a concatenated estimate matrix.

fprintf('Setup CDL channels from Tx to Rx and do perfect channel estimation... \n');

% Normalize rays to always be a cell array where rays{i} contains ray objects for BS i
% Single BS: raytrace returns {[ray objects]} -> already correct
% Two BS: raytrace returns {[ray objects BS1]; [ray objects BS2]} -> already correct
numBS = numel(rays);

channelData = struct();
hEstList = {};


% Loop through each distributed transmitter
for i = 1:numBS

    % rays{i} is always the ray objects for BS i
    % Wrap in cell for createCDLChannel which expects rays{1}
    currentRays = rays(i);
    currentTx = tx(i);

    chFilter = 0;
    %% Try with antennas provided in channel not from ray tracing
    tmpChannel = createCDLChannel(currentRays, fc, currentTx, rx, UEOrientation, chFilter, numRB, SCS);
    % tmpChannel = createCDLChannel2(currentRays, fc, currentTx, rx, UEOrientation, chFilter, numRB, SCS);

    [pathGains, sampleTimes] = tmpChannel();

    pathFilters = getPathFilters(tmpChannel);
    [offset, ~] = nrPerfectTimingEstimate(pathGains, pathFilters);

    % Perfect Channel Estimation
    % hEst dimensions: [Subcarriers, Symbols, RxAnt, TxAnt]
    % Channel tensor of estimated channel across
    % H(f) = sum(a_k * exp(-j*2*pi*f*tau))
    % with a_k = mag and tau = propagation delay
    % (f = subcarriers, t = OFDM symbols (time), Rx-Ant, Tx-Ant)
    % Example: hEst(100, 5, 2, 80) is a single complex double (mag, phase) == specific channel path from Transmit Antenna #80 to Receive Antenna #2, measured at Subcarrier #100 during Time Symbol #5
    % When open in workspace ~ 2D "slice" of the 4D matrix
    % hEst(:,:,1,1) is the 2940x14 matrix representing the channel only between Tx Antenna #1 and Rx Antenna #1, measured across all 2940 subcarriers and 14 time symbols
    % Perfect channel estimate, returned as an N_SC-by-N_SYM-by-N_R-by-N_T complex array, where:
    % N_SC is the number of subcarriers
    % N_SYM is the number of OFDM symbols
    % N_R is the number of receive antennas
    % N_T is the number of transmit antennas
    hEstimation = nrPerfectChannelEstimate(pathGains, pathFilters, numRB, SCS, numSlot, offset, sampleTimes);

    % Use generic naming: BS1, BS2, ... not names
    bsName = "BS" + i;
    channelData.(bsName).Channel = tmpChannel;
    channelData.(bsName).PathGains = pathGains;
    channelData.(bsName).hEst = hEstimation;
    channelData.(bsName).TimingOffset = offset;

    hEstList{end+1} = hEstimation;

    % Plot the mag of the CIR and the timing offset estimate
    % [Nh,Nr] = size(mag);
    % plot(0:(Nh-1),mag,'o:');
    % hold on;
    % plot([offset offset],[0 max(mag(:))*1.25],'k:','LineWidth',2);
    % axis([0 Nh-1 0 max(mag(:))*1.25]);
    % legends = "|h|, antenna " + num2cell(1:Nr);
    % legend([legends "Timing offset estimate"]);
    % ylabel('|h|');
    % xlabel('Channel Impulse Response Samples');
end

% Concatenate for Beamforming
% Stacks matrices along the 4th dimension (Tx Antennas)
% Result: [Subcarriers, Symbols, RxAnt, (TxAnt1 + TxAnt2 + ...)]
% Combine both channels hEst into a "collective Channel" eq 6.1 and 6.5
% p.354, Foundations of User-Centric Cell-Free Massive MIMO, Björnson
% concatenate them (stack them side-by-side)
% hEst is: [Subcarriers, Symbols, Rx, Tx] -> [..., ..., 4, 64]
% cat(4, ...) concatenate other 4 dimensioin -> numElementsTx
hEstCollective = cat(4, hEstList{:});
end
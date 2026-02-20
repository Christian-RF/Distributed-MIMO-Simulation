function [channelData, hEstCollective] = setupDistributedChannels(rays, tx, rx, fc, UEOrientation, numRB, SCS, numSlot, enableDMIMO)
%SETUPDISTRIBUTEDCHANNELS Sets up CDL channels and calculates perfect estimates
%   Returns a struct containing channel objects and a concatenated estimate matrix.

fprintf('Setup CDL channels from Tx to Rx and do perfect channel estimation... \n');

% Setup names and loop targets
if enableDMIMO == 1
    channelNames = ["channelSuperC", "channelAcademica"];
else
    channelNames = "channelSuperC";
end

channelData = struct();
hEstList = {}; % Temporary storage for concatenation

% Loop through each distributed transmitter
for i = 1:numel(channelNames)
    currentName = channelNames(i);

    % Select specific rays and Tx for this link
    if enableDMIMO == 1
        currentRays = rays(i);
        currentTx = tx(i);
    else
        currentRays = rays;
        currentTx = tx;
    end

    % Create Channel Object
    % Channel filter = 0 initially just to get path gains
    chFilter = 0;
    tmpChannel = createCDLChannel(currentRays, fc, currentTx, rx, UEOrientation, chFilter, numRB, SCS, false);

    % Get Path Gains (Snapshot for estimation)
    [pathGains, sampleTimes] = tmpChannel();

    % Perfect Timing Estimation
    pathFilters = getPathFilters(tmpChannel);
    [offset, mag] = nrPerfectTimingEstimate(pathGains, pathFilters);

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

    % Store Data
    channelData.(currentName).Channel = tmpChannel; % Save object for later PDSCH tx
    channelData.(currentName).PathGains = pathGains;
    channelData.(currentName).hEst = hEstimation;
    channelData.(currentName).TimingOffset = offset;

    % Store for concatenation
    hEstList{end+1} = hEstimation;
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
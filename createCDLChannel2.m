function channel = createCDLChannel2(rays, fc, tx, rx, UEOrientation, chFilter, numRB, SCS)

pathToAs = [rays{1}.PropagationDelay] - min([rays{1}.PropagationDelay]);
avgPathGains = -[rays{1}.PathLoss];
pathAoDs = [rays{1}.AngleOfDeparture];
pathAoAs = [rays{1}.AngleOfArrival];
isLOS = any([rays{1}.LineOfSight]);

channel = nrCDLChannel;
channel.DelayProfile = 'Custom';
channel.PathDelays = pathToAs;
channel.AveragePathGains = avgPathGains;
channel.AnglesAoD = pathAoDs(1,:);
channel.AnglesZoD = 90 - pathAoDs(2,:);
channel.AnglesAoA = pathAoAs(1,:);
channel.AnglesZoA = 90 - pathAoAs(2,:);
channel.HasLOSCluster = isLOS;
channel.CarrierFrequency = fc;
channel.NormalizeChannelOutputs = false;
channel.NormalizePathGains = false;

% Channel uses its OWN arrays — not from txsite/rxsite
% TX: 8×8 cross-pol = 128 ports (directional, realistic BS)
channel.TransmitAntennaArray.Size = [8 8 2 1 1];
channel.TransmitAntennaArray.Element = '38.901';
channel.TransmitAntennaArray.ElementSpacing = [0.5 0.5 1 1];
channel.TransmitAntennaArray.PolarizationAngles = [45 -45];
channel.TransmitArrayOrientation = [tx.AntennaAngle(1); -1*tx.AntennaAngle(2); 0];

% RX: 2×1 isotropic cross-pol = 4 ports (omnidirectional, good for dMIMO)
channel.ReceiveAntennaArray.Size = [2 1 2 1 1];
channel.ReceiveAntennaArray.Element = 'isotropic';
channel.ReceiveAntennaArray.ElementSpacing = [0.5 0.5 1 1];
channel.ReceiveAntennaArray.PolarizationAngles = [45 -45];
channel.ReceiveArrayOrientation = [UEOrientation(1); -1*UEOrientation(2); UEOrientation(3)];

ofdmInfo = nrOFDMInfo(numRB, SCS);
channel.SampleRate = ofdmInfo.SampleRate;
channel.ChannelFiltering = chFilter;


end

% % Visualize TX array and cluster paths
% displayChannel(channel, 'LinkEnd', 'Tx');
% 
% % Visualize RX array and cluster paths
% displayChannel(channel, 'LinkEnd', 'Rx');



% %% Verification
% % Verify channel dimensions
% channel.ChannelFiltering = false;
% [pathGains, sampleTimes] = channel();
% fprintf('Path gains size: %s\n', mat2str(size(pathGains)));
% % Expected: [1, numPaths, 4, 128]  (4 RX ports, 128 TX ports)
% 
% % Verify rank from channel estimate
% hEst = nrPerfectChannelEstimate(pathGains, getPathFilters(channel), ...
%     numRB, SCS, 0, 0, sampleTimes);
% fprintf('hEst size: %s\n', mat2str(size(hEst)));
% % Expected: [2940, 14, 4, 128]
% 
% % Check effective rank
% [~, ~, R, P] = size(hEst);
% H = permute(mean(reshape(hEst, [], R, P)), [2 3 1]);
% singVals = svd(H);
% fprintf('Singular values: ');
% fprintf('%.4f ', singVals(1:min(8,end)));
% fprintf('\n');
% fprintf('Effective rank (>1%% of max): %d\n', sum(singVals > 0.01*singVals(1)));
% 
% singVals = svd(H);
% fprintf('Singular values (dB): ');
% fprintf('%.2f ', 20*log10(singVals(1:min(8,end))));
% fprintf('\n');
% fprintf('Condition number: %.2f dB\n', 20*log10(singVals(1)/singVals(end)));
% fprintf('Effective rank (>1%% of max): %d\n', sum(singVals > 0.01*singVals(1)));
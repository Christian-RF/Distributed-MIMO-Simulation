function channel = createCDLChannel(rays, fc, tx, rx, UEOrientation, chFilter, numRB, SCS, false)
%CREATECDLCHANNEL Creates a nrCDchannel object after 3GPP TR 38.901 based on
% raytracing. Based on this example

% (ToA = Time of Arrival)
pathToAs = [rays{1}.PropagationDelay]-min([rays{1}.PropagationDelay]); % Time of arrival of each ray (normalized to 0 sec) 
avgPathGains  = -[rays{1}.PathLoss]; % Average path gains of each ray, -1 because invert pathloss = Gain (dB)

pathAoDs = [rays{1}.AngleOfDeparture]; % Angle of Departure = AoD of each ray (°)
pathAoAs = [rays{1}.AngleOfArrival]; % Angle of Arrival = AoA of each ray (°)

isLOS = any([rays{1}.LineOfSight]); % Line of sight flag

% Set properties of channel
channel = nrCDLChannel;

channel.DelayProfile = 'Custom'; % Set to 'costum' to specify path delays 
channel.PathDelays = pathToAs; % Time of Arrivals
channel.AveragePathGains = avgPathGains; % Path Gains

% Set departure and arrival angles to work from channel (zenith) to rays
% (elevation)
channel.AnglesAoD = pathAoDs(1, :); % azimuth of departure (AoD)
channel.AnglesZoD = 90 - pathAoDs(2, :); % channel uses zenith angle, rays use elevation
channel.AnglesAoA = pathAoAs(1, :); % azimuth of arrival (AoA)
channel.AnglesZoA = 90 - pathAoAs(2, :); % channel uses zenith angle, rays use elevation

channel.HasLOSCluster = isLOS;

channel.CarrierFrequency = fc;

channel.NormalizeChannelOutputs = false; % do not normalize by the number of receive antennas, this would change the receive power
channel.NormalizePathGains = false; % set to false to retain the path gains

% Set channel Tx and Rx antenna arrays
% Convert elevation which is used by the rays to zenith (-1) used by the
% channel -> (-1)
channel.AnglesZoA = 90 - pathAoAs(2, :);    % channel uses zenith angle, rays use elevation
channel.ReceiveAntennaArray = rx.Antenna;
channel.ReceiveArrayOrientation = [UEOrientation(1); -1 * UEOrientation(2); UEOrientation(3)];  % (-1) converts elevation to downtilt

channel.TransmitAntennaArray = tx.Antenna;
channel.TransmitArrayOrientation = [tx.AntennaAngle(1); -1 * tx.AntennaAngle(2); 0];  % (-1) converts elevation to downtilt, but it is already downtilted?


% Set channel sampling rate
ofdmInfo = nrOFDMInfo(numRB, SCS);

channel.SampleRate = ofdmInfo.SampleRate;


% Assmue prefect channel estimation for now (chFilter == 0 == False)
channel.ChannelFiltering = chFilter;

% With no channel filter get channel path gains without sending sig through channel 
% [pathGains,sampleTimes] = channel();


% pg=permute(pathGains,[2 1 3 4]); % first dimension is the number of paths
% if isLOS
%     % in LOS cases sum the first to paths, they correspond to the LOS ray
%     pg = [sum(pg(1:2,:,:,:)); pg(3:end,:,:,:)];
% end
% 
% pg = abs(pg).^2;
% 
% figure;
% plot(pow2db(pg(:,1,1,1)),'o-.');hold on
% plot(avgPathGains,'x-.');hold off
% legend("Instantaneous (1^{st} tx - 1^{st} rx antenna)","Average (from ray tracing)")
% xlabel("Path number"); ylabel("Gain (dB)")
% title('Path gains')

end
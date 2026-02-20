function exposure = calcPower(tx, rx, pm, rays, fc, allocatedTxPower)
tic
fprintf('Starting calculation of coherent and incoherent power/efield ... \n')
exposure = struct();

for i = 1:numel(rays)
    % Get AoA/AoD of each ray
    AoD = [rays{i}.AngleOfDeparture];
    AoA = [rays{i}.AngleOfArrival];

    % Tx -> AoD azimuth and elevation
    azTxGlob = AoD(1,:);
    elTxGlob = AoD(2,:);

    % Change global coordinates (ray tracing) to local (antenna pattern)
    azTx = tx(i).AntennaAngle(1);
    elTx = tx(i).AntennaAngle(2);
    azTxLocal = wrapTo180(azTxGlob - azTx);
    elTxLocal = elTxGlob - elTx;
    elTxLocal = max(min(elTxLocal,90),-90);

    % Rx -> AoA azimuth and elevation
    azRxGlob = AoA(1,:);
    elRxGlob = AoA(2,:);

    % Change global coordinates (ray tracing) to local (antenna pattern)
    % Just one Rx for now
    azRx = rx.AntennaAngle(1);
    elRx = rx.AntennaAngle(2);
    azRxLocal = wrapTo180(azRxGlob - azRx);
    elRxLocal = elRxGlob - elRx;
    elRxLocal = max(min(elRxLocal,90),-90);

    % Calculate gain from AoD(Tx)/AoA(Rx) for each ray
    numRays = length(rays{i});

    gainTx_dBi = zeros(1, numRays);
    gainRx_dBi = zeros(1, numRays);

    for j = 1:numRays
        gainTx_dBi(j) = pattern(tx(i).Antenna, fc, azTxLocal(j), elTxLocal(j));
        gainRx_dBi(j) = pattern(rx.Antenna, fc, azRxLocal(j), elRxLocal(j));
    end

    % Calculate coherent(sigstrength as verification) and incoherent power
    powerTx_dBm = allocatedTxPower;
    pathLoss = [rays{i}.PathLoss];
    phase = [rays{i}.PhaseShift];

    fieldGain_dB = pathLoss - gainTx_dBi - gainRx_dBi;

    % phasor = amplitude * exp(1j * phase);
    amplitude = 10.^(-fieldGain_dB/20);
    phasor = amplitude.*exp(1j.*phase);

    % Calc with isotropic receiver for Efield in dBuv/m
    amplitudeIso = 10.^(-(fieldGain_dB + gainRx_dBi)/20);
    phasorIso = amplitudeIso .* exp(1j .* phase);
    phasorIso_sum = sum(phasorIso);

    exposure.coherentPower_dBm(i) = powerTx_dBm + 20*log10(abs(sum(phasor))); % Voltages -> 20log10

    % Derivation of E-field (dBuV/m) from Isotropic Power (dBm)
    % Power Density S = E^2 / Z0
    % Power Received Prx = S * Aeff
    % Effective Aperture (Isotropic): Aeff = lambda^2 / (4*pi)
    % Prx = (E^2 / Z0) * (lambda^2 / (4*pi))
    % E = sqrt( Prx * Z0 * 4*pi ) / lambda
    % Convert to dB: 20*log10(E) = 10*log10(Prx) + 10*log10(Z0*4*pi) - 20*log10(lambda)
    %
    % Prx is in dBm, subtract 30 to get dBW (Watts)
    % E is in dBuV, add 120 to convert from dBV (Volts)
    exposure.choherentEfield_dbuv(i) = powerTx_dBm...
        + 20*log10(abs(phasorIso_sum)) - 30 ...
        + 10*log10(rfprop.Constants.Z0 * 4 * pi)...
        + 120 ...
        - 20*log10(rfprop.Constants.LightSpeed / fc);

    exposure.incoherentPower_dBm(i) = powerTx_dBm + 10*log10(sum(abs(phasor).^2)); % Power -> 10log10
    % exposure.incoherentPower_dBm = powerTx_dBm + 10*log10(sum(abs(amplitude).^2));
    % exposure.sigstrengthPower_dBm(i) = sigstrength(rx, tx(i), pm, Type = 'power');
    % exposure.sigstrengthEfield_dBuv(i) = sigstrength(rx, tx(i), pm, Type = 'efield');
end

% Sigstrength is vectorized for multiple tx but its not much faster then doing
% it in the loop
% In loop: Elapsed time is 108.186635 seconds
% Out of loop: Elapsed time is 104.345622 seconds.
% exposure.sigstrengthPower_dBm = sigstrength(rx, tx, pm, Type = 'power')';
% exposure.sigstrengthEfield_dBuv = sigstrength(rx, tx, pm, Type = 'efield')';

fprintf('With a duration of \n');
toc
% % % % Get AoA/AoD of each ray
% % % AoD = [rays{i}.AngleOfDeparture];
% % % AoA = [rays{1}.AngleOfArrival];
% % % 
% % % % Tx -> AoD azimuth and elevation
% % % azTxGlob = AoD(1,:);
% % % elTxGlob = AoD(2,:);
% % % 
% % % % Change global coordinates (ray tracing) to local (antenna pattern)
% % % azTx = tx.AntennaAngle(1);
% % % elTx = tx.AntennaAngle(2);
% % % azTxLocal = wrapTo180(azTxGlob - azTx);
% % % elTxLocal = elTxGlob - elTx;
% % % elTxLocal = max(min(elTxLocal,90),-90);
% % % 
% % % % Rx -> AoA azimuth and elevation
% % % azRxGlob = AoA(1,:);
% % % elRxGlob = AoA(2,:);
% % % 
% % % % Change global coordinates (ray tracing) to local (antenna pattern)
% % % azRx = rx.AntennaAngle(1);
% % % elRx = rx.AntennaAngle(2);
% % % azRxLocal = wrapTo180(azRxGlob - azRx);
% % % elRxLocal = elRxGlob - elRx;
% % % elRxLocal = max(min(elRxLocal,90),-90);
% % % 
% % % % Calculate gain from AoD(Tx)/AoA(Rx) for each ray
% % % numRays = length(rays{1});
% % % 
% % % gainTx_dBi = zeros(1, numRays);
% % % gainRx_dBi = zeros(1, numRays);
% % % 
% % % for i = 1:numRays
% % %     gainTx_dBi(i) = pattern(tx.Antenna, fc, azTxLocal(i), elTxLocal(i));
% % %     gainRx_dBi(i) = pattern(rx.Antenna, fc, azRxLocal(i), elRxLocal(i));
% % % end
% % % 
% % % % Calculate coherent(sigstrength as verification) and incoherent power
% % % powerTx_dBm = allocatedTxPower;
% % % pathLoss = [rays{1}.PathLoss];
% % % phase = [rays{1}.PhaseShift];
% % % 
% % % fieldGain_dB = pathLoss - gainTx_dBi - gainRx_dBi;
% % % 
% % % % phasor = amplitude * exp(1j * phase);
% % % amplitude = 10.^(-fieldGain_dB/20);
% % % phasor = amplitude.*exp(1j.*phase);
% % % 
% % % % Calc with isotropic receiver for Efield in dBuv/m
% % % amplitudeIso = 10.^(-(fieldGain_dB + gainRx_dBi)/20);
% % % phasorIso = amplitudeIso .* exp(1j .* phase);
% % % phasorIso_sum = sum(phasorIso);
% % % 
% % % exposure.coherentPower_dBm = powerTx_dBm + 20*log10(abs(sum(phasor))); % Voltages -> 20log10
% % % 
% % % % Derivation of E-field (dBuV/m) from Isotropic Power (dBm)
% % % % Power Density S = E^2 / Z0
% % % % Power Received Prx = S * Aeff
% % % % Effective Aperture (Isotropic): Aeff = lambda^2 / (4*pi)
% % % % Prx = (E^2 / Z0) * (lambda^2 / (4*pi))
% % % % E = sqrt( Prx * Z0 * 4*pi ) / lambda
% % % % Convert to dB: 20*log10(E) = 10*log10(Prx) + 10*log10(Z0*4*pi) - 20*log10(lambda)
% % % %
% % % % Prx is in dBm, subtract 30 to get dBW (Watts)
% % % % E is in dBuV, add 120 to convert from dBV (Volts)
% % % exposure.choherentEfield_dbuv = powerTx_dBm...
% % %     + 20*log10(abs(phasorIso_sum)) - 30 ...
% % %     + 10*log10(rfprop.Constants.Z0 * 4 * pi)...
% % %     + 120 ...
% % %     - 20*log10(rfprop.Constants.LightSpeed / fc);
% % % 
% % % exposure.incoherentPower_dBm = powerTx_dBm + 10*log10(sum(abs(phasor).^2)); % Power -> 10log10
% % % % exposure.incoherentPower_dBm = powerTx_dBm + 10*log10(sum(abs(amplitude).^2));
% % % exposure.sigstrengthPower_dBm = sigstrength(rx, tx, pm, Type = 'power');
% % % exposure.sigstrengthEfield_dBuv = sigstrength(rx, tx, pm, Type = 'efield');
% % % 

function exposure = calcPower4(tx, rx, rays, fc, allocatedTxPower, gainTx, gainRx)

% Preallocate
numBS = numel(rays);
exposure = struct();
exposure.coherentPower_dBm = -inf(1, numBS);
exposure.incoherentPower_dBm = -inf(1, numBS);
exposure.coherentPowerIso_dBm = -inf(1, numBS);
exposure.incoherentPowerIso_dBm = -inf(1, numBS);
exposure.coherentEfield_dBuv = -inf(1, numBS);
exposure.incoherentEfield_dBuv = -inf(1, numBS);
allPhasors = [];
allPhasorsIso = [];

for i = 1:numBS
    AoD = [rays{i}.AngleOfDeparture];
    AoA = [rays{i}.AngleOfArrival];

    % Tx local angles
    azTxLocal = wrapTo180(AoD(1,:) - tx(i).AntennaAngle(1));
    elTxLocal = AoD(2,:) - tx(i).AntennaAngle(2);
    elTxLocal = max(min(elTxLocal, 90), -90);

    % Rx local angles
    azRxLocal = wrapTo180(AoA(1,:) - rx.AntennaAngle(1));
    elRx = 0;
    if numel(rx.AntennaAngle) >= 2
        elRx = rx.AntennaAngle(2);
    end
    elRxLocal = AoA(2,:) - elRx;
    elRxLocal = max(min(elRxLocal, 90), -90);

    % interp2 instead of pattern()
    % for j = 1:numRays
    %     gainTx_dBi(j) = pattern(tx(i).Antenna, fc, azTxLocal(j), elTxLocal(j));
    %     %gainRx_dBi(j) = pattern(rx.Antenna, fc, azRxLocal(j), elRxLocal(j));
    %
    %     % Rx antenna: handle isotropic (string) vs actual antenna object
    %     if isstring(rx.Antenna) || ischar(rx.Antenna)
    %         gainRx_dBi(j) = 0;  % isotropic = 0 dBi everywhere for sampling exposure map
    %     else
    %         gainRx_dBi(j) = pattern(rx.Antenna, fc, azRxLocal(j), elRxLocal(j));
    %     end
    % end
    % gain(i) picks the correct BS gain grid
    gainTx_dBi = interp2(gainTx(i).az, gainTx(i).el, gainTx(i).Grid, azTxLocal, elTxLocal, 'linear');

    if isstring(rx.Antenna) || ischar(rx.Antenna)
        gainRx_dBi = zeros(1, length(rays{i}));
    else
        gainRx_dBi = interp2(gainRx.az, gainRx.el, gainRx.rxGrid, azRxLocal, elRxLocal, 'linear');
    end

    powerTx_dBm = allocatedTxPower(i); % Power per basestation
    pathLoss = [rays{i}.PathLoss];
    phase = [rays{i}.PhaseShift];

    fieldGain_dB = pathLoss - gainTx_dBi - gainRx_dBi;

    amplitude = 10.^(-fieldGain_dB/20);
    phasor = amplitude .* exp(1j .* phase);

    % Isotropic phasor (no Rx gain) for E-field
    amplitudeIso = 10.^(-(fieldGain_dB + gainRx_dBi)/20);
    phasorIso = amplitudeIso .* exp(1j .* phase);
    phasorIso_sum = sum(phasorIso);

    % Needed for different resource allocations
    sqrtPwr = sqrt(10^((powerTx_dBm - 30) / 10));
    allPhasors    = [allPhasors,    phasor * sqrtPwr];
    allPhasorsIso = [allPhasorsIso, phasorIso * sqrtPwr];

    % Power with Rx gain (for SNR, link budget)
    exposure.coherentPower_dBm(i) = powerTx_dBm + 20*log10(abs(sum(phasor)));
    exposure.incoherentPower_dBm(i) = powerTx_dBm + 10*log10(sum(abs(phasor).^2));

    % Power without Rx gain (isotropic, for E-field derivation)
    exposure.coherentPowerIso_dBm(i) = powerTx_dBm + 20*log10(abs(phasorIso_sum));
    exposure.incoherentPowerIso_dBm(i) = powerTx_dBm + 10*log10(sum(abs(phasorIso).^2));

    % Per-BS E-field (from isotropic power)
    exposure.coherentEfield_dBuv(i) = exposure.coherentPowerIso_dBm(i) - 30 ...
        + 10*log10(rfprop.Constants.Z0 * 4 * pi) ...
        + 120 ...
        - 20*log10(rfprop.Constants.LightSpeed / fc);

    exposure.incoherentEfield_dBuv(i) = exposure.incoherentPowerIso_dBm(i) - 30 ...
        + 10*log10(rfprop.Constants.Z0 * 4 * pi) ...
        + 120 ...
        - 20*log10(rfprop.Constants.LightSpeed / fc);
end

%% Global coherent
exposure.totalCoherentPower_dBm   = 10*log10(abs(sum(allPhasors))^2) + 30;
exposure.totalIncoherentPower_dBm = 10*log10(sum(abs(allPhasors).^2)) + 30;


%% Global E-field from isotropic phasors
exposure.totalCoherentEfield_dBuv = 20*log10(abs(sum(allPhasorsIso))) ...
    + 10*log10(rfprop.Constants.Z0 * 4 * pi) ...
    + 120 ...
    - 20*log10(rfprop.Constants.LightSpeed / fc);

exposure.totalIncoherentEfield_dBuv = 10*log10(sum(abs(allPhasorsIso).^2)) ...
    + 10*log10(rfprop.Constants.Z0 * 4 * pi) ...
    + 120 ...
    - 20*log10(rfprop.Constants.LightSpeed / fc);
end
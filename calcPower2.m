function exposure = calcPower(tx, rx, rays, fc, allocatedTxPower)
tic
fprintf('Starting calculation of coherent and incoherent power/efield ... \n')
exposure = struct();

for i = 1:numel(rays)
    AoD = [rays{i}.AngleOfDeparture];
    AoA = [rays{i}.AngleOfArrival];

    azTxGlob = AoD(1,:);
    elTxGlob = AoD(2,:);
    azTx = tx(i).AntennaAngle(1);
    elTx = tx(i).AntennaAngle(2);
    azTxLocal = wrapTo180(azTxGlob - azTx);
    elTxLocal = elTxGlob - elTx;
    elTxLocal = max(min(elTxLocal,90),-90);

    azRxGlob = AoA(1,:);
    elRxGlob = AoA(2,:);
    azRx = rx.AntennaAngle(1);
    elRx = rx.AntennaAngle(2);
    azRxLocal = wrapTo180(azRxGlob - azRx);
    elRxLocal = elRxGlob - elRx;
    elRxLocal = max(min(elRxLocal,90),-90);

    numRays = length(rays{i});
    gainTx_dBi = zeros(1, numRays);
    gainRx_dBi = zeros(1, numRays);

    for j = 1:numRays
        gainTx_dBi(j) = pattern(tx(i).Antenna, fc, azTxLocal(j), elTxLocal(j), ...
            'Type', 'powerdb', 'Normalize', false);
        gainRx_dBi(j) = pattern(rx.Antenna, fc, azRxLocal(j), elRxLocal(j), ...
            'Type', 'powerdb', 'Normalize', false);
    end

    powerTx_dBm = allocatedTxPower;
    pathLoss = [rays{i}.PathLoss];
    phase = [rays{i}.PhaseShift];

    fieldGain_dB = pathLoss - gainTx_dBi - gainRx_dBi;

    amplitude = 10.^(-fieldGain_dB/20);
    phasor = amplitude .* exp(1j .* phase);

    amplitudeIso = 10.^(-(fieldGain_dB + gainRx_dBi)/20);
    phasorIso = amplitudeIso .* exp(1j .* phase);
    phasorIso_sum = sum(phasorIso);

    exposure.coherentPower_dBm(i) = powerTx_dBm + 20*log10(abs(sum(phasor)));
    exposure.incoherentPower_dBm(i) = powerTx_dBm + 10*log10(sum(abs(phasor).^2));

    exposure.coherentPowerIso_dBm(i) = powerTx_dBm + 20*log10(abs(phasorIso_sum));
    exposure.incoherentPowerIso_dBm(i) = powerTx_dBm + 10*log10(sum(abs(phasorIso).^2));

    exposure.coherentEfield_dBuv(i) = exposure.coherentPowerIso_dBm(i) - 30 ...
        + 10*log10(rfprop.Constants.Z0 * 4 * pi) ...
        + 120 ...
        - 20*log10(rfprop.Constants.LightSpeed / fc);

    exposure.incoherentEfield_dBuv(i) = exposure.incoherentPowerIso_dBm(i) - 30 ...
        + 10*log10(rfprop.Constants.Z0 * 4 * pi) ...
        + 120 ...
        - 20*log10(rfprop.Constants.LightSpeed / fc);
end

exposure.totalIncoherentPower_dBm = 10*log10(sum(10.^(exposure.incoherentPower_dBm/10)));
exposure.totalCoherentPower_dBm = 10*log10(sum(10.^(exposure.coherentPower_dBm/10)));

exposure.totalIncoherentEfield_dBuv = 20*log10(sqrt(sum(10.^(exposure.incoherentPowerIso_dBm/10)))) ...
    - 30 + 10*log10(rfprop.Constants.Z0 * 4 * pi) + 120 ...
    - 20*log10(rfprop.Constants.LightSpeed / fc);
exposure.totalCoherentEfield_dBuv = 20*log10(sqrt(sum(10.^(exposure.coherentPowerIso_dBm/10)))) ...
    - 30 + 10*log10(rfprop.Constants.Z0 * 4 * pi) + 120 ...
    - 20*log10(rfprop.Constants.LightSpeed / fc);

fprintf('With a duration of \n');
toc
end
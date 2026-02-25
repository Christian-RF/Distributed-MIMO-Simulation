function exposure = calcPowerPol(tx, rx, rays, fc, allocatedTxPower, xpr_dB)
% Full polarimetric power calculation
% xpr_dB: XPR from TR 38.901 Table 7.5-6 (e.g., 8 for UMi NLOS)

    N_rays = numel(rays);
    exposure = struct();
    exposure.coherentPower_dBm      = -inf(1, N_rays);
    exposure.incoherentPower_dBm    = -inf(1, N_rays);
    exposure.coherentPowerIso_dBm   = -inf(1, N_rays);
    exposure.incoherentPowerIso_dBm = -inf(1, N_rays);
    
    if nargin < 6
        xpr_dB = 8;  % default: UMi NLOS
    end

    for i = 1:numel(rays)
        AoD = [rays{i}.AngleOfDeparture];
        AoA = [rays{i}.AngleOfArrival];
        
        % --- Local angle transforms (same as before) ---
        azTx = tx(i).AntennaAngle(1);
        elTx = tx(i).AntennaAngle(2);
        azTxLocal = wrapTo180(AoD(1,:) - azTx);
        elTxLocal = max(min(AoD(2,:) - elTx, 90), -90);
        
        azRx = rx.AntennaAngle(1);
        elRx = 0;
        if numel(rx.AntennaAngle) >= 2
            elRx = rx.AntennaAngle(2);
        end
        azRxLocal = wrapTo180(AoA(1,:) - azRx);
        elRxLocal = max(min(AoA(2,:) - elRx, 90), -90);
        
        numRays = length(rays{i});
        
        % --- Get polarized field patterns ---
        [Ftx_theta, Ftx_phi] = getPolComponents(tx(i).Antenna, fc, azTxLocal, elTxLocal);
        [Frx_theta, Frx_phi] = getPolComponents(rx.Antenna,   fc, azRxLocal, elRxLocal);
        
        % --- Path parameters ---
        pathLoss  = [rays{i}.PathLoss];        % dB
        phase     = [rays{i}.PhaseShift];       % rad
        amplitude = 10.^(-pathLoss/20);         % linear voltage
        
        % --- Generate polarization matrices ---
        P = genPolMatrix(xpr_dB, numRays);
        
        % --- Compute polarimetric phasor per ray ---
        phasor    = zeros(1, numRays);
        phasorIso = zeros(1, numRays);
        
        for j = 1:numRays
            % Tx field vector [2x1]
            f_tx = [Ftx_theta(j); Ftx_phi(j)];
            
            % Rx field vector [2x1]
            f_rx = [Frx_theta(j); Frx_phi(j)];
            
            % Polarimetric channel coefficient for this ray
            %   h = f_rx^T * P * f_tx * amplitude * exp(j*phase)
            h_pol = (f_rx.' * P{j} * f_tx);
            
            phasor(j) = h_pol * amplitude(j) * exp(1j * phase(j));
            
            % Isotropic Rx: replace f_rx with unit response
            %   Isotropic in both pols → [1/sqrt(2); 1/sqrt(2)]
            f_rx_iso = [1; 1] / sqrt(2);
            h_pol_iso = (f_rx_iso.' * P{j} * f_tx);
            phasorIso(j) = h_pol_iso * amplitude(j) * exp(1j * phase(j));
        end
        
        % --- Power calculations ---
        powerTx_dBm = allocatedTxPower;
        
        exposure.coherentPower_dBm(i)   = powerTx_dBm + 20*log10(abs(sum(phasor)));
        exposure.incoherentPower_dBm(i) = powerTx_dBm + 10*log10(sum(abs(phasor).^2));
        
        exposure.coherentPowerIso_dBm(i)   = powerTx_dBm + 20*log10(abs(sum(phasorIso)));
        exposure.incoherentPowerIso_dBm(i) = powerTx_dBm + 10*log10(sum(abs(phasorIso).^2));
    end
    
    % --- Totals ---
    exposure.totalIncoherentPower_dBm = 10*log10(sum(10.^(exposure.incoherentPower_dBm/10)));
    exposure.totalCoherentPower_dBm   = 10*log10(sum(10.^(exposure.coherentPower_dBm/10)));
end
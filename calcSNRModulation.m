function [modOrder, usedMod, SNR_dB, coderate, SE] = calcSNRModulation(exposureAfterBF, allocBW, SCS, MCSIndexTable)
%%CALCSNRMODULATION calculates SNR, SE, Modulation and coderate based MCS table 
% and the provided signal strength in exposure

%% Calc thermal noise over BW (Subcarriers)
k = physconst('Boltzmann'); % 1.38e-23 (J/K)
T = 290; % Temp in Kelvin
rxNF = 9; % RX-Noise-Figure (dB)
% n = k*T*BW
% SCS in kHz -> *1e3 
% n0 = k * T * SCS * 1e3 * 10^(rxNF / 10); % Thermal Noise in W of one subcarrier 30kHz
% Not per single carrier but noise in complete band -> scs*2940 = allocBW
% n0 = k * T * SCS * 1e3 * (allocBW/(SCS*1e3)) * 10^(rxNF / 10);
n0 = k * T * allocBW * 10^(rxNF / 10);
n0_mW = n0 * 1e3; % Convert to mW for correct calculation with sigStrength

%% Calc avarage Rx power per Subcarrier
% Not per subcarrier anymore since incoherent power approximates power in
% complete band
sigStrength_mW = 10.^(exposureAfterBF.incoherentPower_dBm/10);% dBm to mW
% avgSigStrength = sum(sigStrength_mW) / length(sigStrength_mW);  % Average receive Power over Subcarrier

%% Calculate SNR and SE
SNR_Lin = sigStrength_mW ./ n0_mW; % should this be n0 or n0_mW?????
% SNR for PDSCH simulation needed -> is very high since there is no
% interference implemented -> maybe set lower manual
SNR_dB = 10.*log10(SNR_Lin);

%% Calc spectral efficiency for finding modulation and coderate
% Find good factor combinations of both and set manual
SE = log2(1 + SNR_Lin); % (bit/(sec Hz))

%% Calc Modulation and Coderate based on SNR and modulation table
numTx = length(SE);

% Pre-allocate arrays for speed
modOrder = zeros(1, numTx);
coderate = zeros(1, numTx);
usedMod = cell(1, numTx); % Use cell array for strings

for i = 1:numTx
    % Find the matching index in the MCS table for the i-th transmitter
    idx = find(MCSIndexTable.Spectral_efficiency >= SE(i), 1);
    
    % If SE is higher than the table max, use the highest available MCS
    if isempty(idx)
        idx = length(MCSIndexTable.Spectral_efficiency);
    end
    
    % Extract values
    modOrder(i) = MCSIndexTable.Modulation_Order(idx);
    coderate(i) = MCSIndexTable.Target_Code_Rate(idx);
    
    % Determine Modulation String
    switch modOrder(i)
        case 2
            usedMod{i} = 'QPSK';
        case 4
            usedMod{i} = '16QAM';
        case 6
            usedMod{i} = '64QAM';
        case 8
            usedMod{i} = '256QAM';
    end
end

end
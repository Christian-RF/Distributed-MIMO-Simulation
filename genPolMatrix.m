function P = genPolMatrix(xpr_dB, numRays)
% Generate 2x2 polarization transfer matrices per TR 38.901
% xpr_dB: cross-polarization ratio in dB (e.g., 8 for urban NLOS)
%         Table 7.5-6 in TR 38.901:
%         UMi-LOS: 9 dB, UMi-NLOS: 8 dB
%         UMa-LOS: 8 dB, UMa-NLOS: 7 dB
%
% Returns: P{j} = 2x2 complex matrix for each ray j

    kappa_inv = 10^(-xpr_dB/10);  % linear, inverse XPR
    sqrtK = sqrt(kappa_inv);
    
    P = cell(1, numRays);
    for j = 1:numRays
        % 4 independent uniform random phases
        phi_tt = 2*pi*rand();
        phi_tp = 2*pi*rand();
        phi_pt = 2*pi*rand();
        phi_pp = 2*pi*rand();
        
        P{j} = [exp(1j*phi_tt),       sqrtK * exp(1j*phi_tp); ...
                sqrtK * exp(1j*phi_pt), exp(1j*phi_pp)];
    end
end
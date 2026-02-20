%% 38.214 - Table 5.1.3.1-2: MCS index table 2 for PDSCH
% QPSK up to 256 QAM

MCSIndexTable = table('Size',[28 5],'VariableTypes',["double", "double","string", "double", "double"],...
                                    'VariableNames',["MCS_Index", "Modulation_Order", "Modulation","Target_Code_Rate", "Spectral_efficiency"]);

MCSIndexTable.MCS_Index = (0:27)';

MCSIndexTable.Modulation_Order(1:5) = 2; % QPSK
MCSIndexTable.Modulation_Order(6:11) = 4; % 16 QAM
MCSIndexTable.Modulation_Order(12:20) = 6; % 64 QAM
MCSIndexTable.Modulation_Order(21:28) = 8; % 256 QAM

MCSIndexTable.Modulation(1:5) = 'QPSK'; % QPSK
MCSIndexTable.Modulation(6:11) = '16QAM'; % 16 QAM
MCSIndexTable.Modulation(12:20) = '64QAM'; % 64 QAM
MCSIndexTable.Modulation(21:28) = '256QAM'; % 256 QAM


MCSIndexTable.Target_Code_Rate(1) = 120;
MCSIndexTable.Target_Code_Rate(2) = 193;
MCSIndexTable.Target_Code_Rate(3) = 308;
MCSIndexTable.Target_Code_Rate(4) = 449;
MCSIndexTable.Target_Code_Rate(5) = 602;
MCSIndexTable.Target_Code_Rate(6) = 378;
MCSIndexTable.Target_Code_Rate(7) = 434;
MCSIndexTable.Target_Code_Rate(8) = 490;
MCSIndexTable.Target_Code_Rate(9) = 553;
MCSIndexTable.Target_Code_Rate(10) = 616;
MCSIndexTable.Target_Code_Rate(11) = 658;
MCSIndexTable.Target_Code_Rate(12) = 466;
MCSIndexTable.Target_Code_Rate(13) = 517;
MCSIndexTable.Target_Code_Rate(14) = 567;
MCSIndexTable.Target_Code_Rate(15) = 616;
MCSIndexTable.Target_Code_Rate(16) = 666;
MCSIndexTable.Target_Code_Rate(17) = 719;
MCSIndexTable.Target_Code_Rate(18) = 772;
MCSIndexTable.Target_Code_Rate(19) = 822;
MCSIndexTable.Target_Code_Rate(20) = 873;
MCSIndexTable.Target_Code_Rate(21) = 682.5;
MCSIndexTable.Target_Code_Rate(22) = 711;
MCSIndexTable.Target_Code_Rate(23) = 754;
MCSIndexTable.Target_Code_Rate(24) = 797;
MCSIndexTable.Target_Code_Rate(25) = 841;
MCSIndexTable.Target_Code_Rate(26) = 885;
MCSIndexTable.Target_Code_Rate(27) = 916.5;
MCSIndexTable.Target_Code_Rate(28) = 948;

MCSIndexTable.Spectral_efficiency(1) = 0.2344;
MCSIndexTable.Spectral_efficiency(2) = 0.3770;
MCSIndexTable.Spectral_efficiency(3) = 0.6016;
MCSIndexTable.Spectral_efficiency(4) = 0.8770;
MCSIndexTable.Spectral_efficiency(5) = 1.1758;
MCSIndexTable.Spectral_efficiency(6) = 1.4766;
MCSIndexTable.Spectral_efficiency(7) = 1.6953;
MCSIndexTable.Spectral_efficiency(8) = 1.9141;
MCSIndexTable.Spectral_efficiency(9) = 2.1602;
MCSIndexTable.Spectral_efficiency(10) = 2.4063;
MCSIndexTable.Spectral_efficiency(11) = 2.5703;
MCSIndexTable.Spectral_efficiency(12) = 2.7305;
MCSIndexTable.Spectral_efficiency(13) = 3.0293;
MCSIndexTable.Spectral_efficiency(14) = 3.3223;
MCSIndexTable.Spectral_efficiency(15) = 3.6094;
MCSIndexTable.Spectral_efficiency(16) = 3.9023;
MCSIndexTable.Spectral_efficiency(17) = 4.2129;
MCSIndexTable.Spectral_efficiency(18) = 4.5234;
MCSIndexTable.Spectral_efficiency(19) = 4.8164;
MCSIndexTable.Spectral_efficiency(20) = 5.1152;
MCSIndexTable.Spectral_efficiency(21) = 5.3320;
MCSIndexTable.Spectral_efficiency(22) = 5.5547;
MCSIndexTable.Spectral_efficiency(23) = 5.8906;
MCSIndexTable.Spectral_efficiency(24) = 6.2266;
MCSIndexTable.Spectral_efficiency(25) = 6.5703;
MCSIndexTable.Spectral_efficiency(26) = 6.9141;
MCSIndexTable.Spectral_efficiency(27) = 7.1602;
MCSIndexTable.Spectral_efficiency(28) = 7.4063;

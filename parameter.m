% PARAMETER provides all needed parameter for the simulation.

%% Base Parameterss
c = physconst('lightspeed'); % (m/s)
fc = 3.655e9; % Center frequency (1/s)
lambda = c/fc; % Wavelength (m)
bw = 90e6; % 90MHz bandwidth
% frequencies = [fc]; % Use fc and bw start/end for evaluation fc-(bw/2) fc fc+(bw/2) take band beginning and end later into acount
thermalNoise = -174; % (dBm/Hz)
EIRP = 44; % (dBm) -> 25W,  40dBm=10W, ITU-R M.2412 Table 5b) 44dBm?

SLA = 30; % (dB) Maximum side-lobe level attenuation

txNF = 5; % (dB) noise figure
rxNF = 7; % (dB) noise figure

txElementGain = 8; % BS antenna element gain (dBi)
rxElementGain = 0; % UE antenna element gain (dBi)

UEOrientation = [0 0 0]; % Orientation of antenna elements


% Define the Antenna Elements (TR 38.901 Table 7.3-1)
% Create two elements to represent the Cross-Pol pair (P=2)
% Both share the same radiation pattern parameters (Gain, Beamwidth, etc.)

elementParameter = { ...
    'Beamwidth', [65 65], ...       % theta_3dB, phi_3dB
    'MaximumGain', 8, ...           % GE,max (8 dBi)
    'SidelobeLevel', 30, ...        % SLAv (30 dB)
    'MaximumAttenuation', 30, ...   % A_max (30 dB)
    'PolarizationModel', 2,...      % Model-2 is standard for TR 38.901
    'FrequencyRange' [3.3e9 3.8e9]};% band n78 or just channel fc-bw/2 fc+bw/2       

%% Comeback later for polarization implementation
% % % Element 1: +45 degree slant
elPlus45 = phased.NRAntennaElement(elementParameter{:}, ...
    'FrequencyRange', [3.3e9 3.8e9],...
    PolarizationAngle = 45);

% Element 2: -45 degree slant
elMinus45 = phased.NRAntennaElement(elementParameter{:}, ...
    'FrequencyRange', [3.3e9 3.8e9],...
    PolarizationAngle = -45);

element = phased.NRAntennaElement(elementParameter{:});

% nrElement = phased.NRAntennaElement('FrequencyRange', [3.3e9 3.8e9]);
% txArrayTest = phased.NRRectangularPanelArray(ElementSet=default);

% Define the Rectangular Panel Array
% Mapping the tuple (Mg, Ng, M, N, P) to MATLAB properties:
% Size = [M, N, Mg, Ng]  <-- Note the order!
% Spacing = [d_v, d_h, dg_v, dg_h]

% Elements per panel (Rows, Cols)
M = 8;
N = 8;

% Number of panels (Vertical, Horizontal)
Mg = 1;
Ng = 1;

% Similar to the default NRRect
% txArray = phased.NRRectangularPanelArray( ... 
%     'ElementSet', {elMinus45, elPlus45}, ...   % Assign the cross-pol pair
%     'Size', [M, N, Mg, Ng], ...                  % 8x8 Array, 1 Panel
%     'Spacing', [0.5, 0.5, 1, 1] * lambda'); % standard spacing 
    % ,...
    % 'EnablePanelSubarray',true,...
    % 'SubarraySteering', 'Custom);

% % Visualization to Verify
% figure;
% viewArray(txArray, 'ShowIndex', 'All');
% title('3GPP TR 38.901 Cross-Polarized Panel Array (8x8)');

% txArray = phased.NRRectangularPanelArray( ... 
%     'ElementSet', {elMinus45, elPlus45}, ...   % Assign the cross-pol pair
%     'Size', [M, N, Mg, Ng], ...                  % 8x8 Array, 1 Panel
%     'Spacing', [0.5, 0.5, 1, 1] * lambda); % standard spacing 

txArray = phased.NRRectangularPanelArray( ... 
    'ElementSet', {element}, ... % Single element for ray tracing
    'Size', [M, N, Mg, Ng], ...  % 8x8 Array, 1 Panel
    'Spacing', [0.5, 0.5, 1, 1] * lambda);

%% Rx Antenna parameter

rxArraySize = [2 2]; % Size of UE array with 2x2 elements, typical smartphone with 4 antennas [2 2]

ueIsotropic = phased.IsotropicAntennaElement(FrequencyRange = [3.3e9 3.8e9]); % Creates a Isotrpic antenna for n78 band

% ueIsotropic = phased.CustomAntennaElement(FrequencyVector = [3.3e9 3.8e9]);
% ueIsotropic.SpecifyPolarizationPattern = true;
rxArray = phased.URA(rxArraySize, ...
    'Element',ueIsotropic,...
    'ElementSpacing',[0.5 0.5] * lambda);

% rxArraySize = [1 2];
% ueElement = phased.CrossedDipoleAntennaElement(...
%     'FrequencyRange', [3.3e9 3.8e9]);
% rxArray = phased.ULA('NumElements', 2, ...
%     'Element',ueElement,...
%     'ElementSpacing',0.5 * lambda);




% figure;
% viewArray(rxArray, 'ShowIndex', 'All');


%% Tx Setup

% First Tx site on the WiWi building (Sammelbau Templergraben 64) near Super C
% gNB ID 83008293 (Macro) - NR
% Standortbescheinigungs-Nr.: 520370
% check again on cellmapper, coordinates could be off and hight is
% important, also another BS is on that roof
txLatWiwi = 50.7776144; % (°)
txLonWiwi = 6.0791105; % (°)

% Second Tx site on mensa academica
% gNB ID 83035161 (Macro) - NR
% Standortbescheinigungs-Nr.: 520142
txLatAcademica = 50.7813726; % (°)
txLonAcademica = 6.0763231; % (°)

% Each bs is 10 meter heigh
bsHeight = 10;


%% Rx Setup

% Rx midpoint between both bs
% rxLat = 50.7794935;
% rxLong = 6.0777168;

% Rx behind superC
% rxLat = 50.7784224;
% rxLong = 6.0783278;

rxLat = 50.7782852; % (°) NLOS 50.7791588 LOS 50.7782852
rxLong = 6.0791742; % (°) NLOS 6.0781231 LOS 6.0791742

% 1.5m height as normal Rx height
ueHeight = 1.5;


%% 5G Parameters
% Traffic model Full Buffer? Max exposition with burst of full? Worstcase

% Bandwidth configuration, required to set the channel sampling rate and for perfect channel estimation
% From TR 38.101 Table 5.3.2-1 max transmission bw for FR1
SCS = 30; % subcarrier spacing
numRB = 245; % number of resource blocks for 90 MHz bandwidth, takes Guardbands into account
numSubcarrier = numRB * 12; % 2940 Subcarriers
guardBand = 885e3; % TS 38.521-1 T5.3.3.-1 (BW*numRB*SCS*12)/2 -SCS/2

numLayers = 4; % Number of data streams (layers) later upto 4?
scOffset = 0; % Subcarrier offset index -> use all subcarriers

numFrames = 5; % Number of Frames/Slots  10
numHARQ = 16; % Number of HARQ = Hybrid Automatic Repeat Request 16
numSlot = 0;


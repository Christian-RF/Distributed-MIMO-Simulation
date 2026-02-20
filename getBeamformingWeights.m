function [weightsTx, weightsRx, D] = getBeamformingWeights(hEstimation, nLayers, scOffset, numRB)
% Get beamforming weights given a channel matrix hEst and the number of
% layers nLayers. One set of weights is provided for the whole bandwidth.
% The beamforming weights are calculated using singular value (SVD)
% decomposition.
%
% Only part of the channel estimate is used to get the weights, this is
% indicated by an offset SCOFFSET (offset from the first subcarrier) and a
% width in RBs (NORBS).

% Average channel estimate
[~, ~ ,R , P] = size(hEstimation); % [~] just placeholder, R and P = numbers of Rx and Tx antennas

scNo = scOffset + 1; % Starting subcarrier idx

hEstimation = hEstimation(scNo:scNo + (12 * numRB-1), :, :, :); % (1:2940, :, :, :) select all subcarriers

% Reshape 4D hEst matrix into 3D matrix
% Calc average of channel values over all subcarriers and timeslots for
% each Rx-Tx antenna pair
% Permute rearranges the matrix
H = permute(mean(reshape(hEstimation, [], R, P)), [2 3 1]); % Rx-Tx 4x128 Antenna

%% SVD decomposition of H
% H = 4x128 Channel matrix

% U = 4x4 (Rx Rx) matrix of left singular vectors, represent the optimal
% "combining" vectors (columns) for Rx

% D = 4x128 matrix of singular values (4) on its diagonal, they represent
% the channel gains of the independet spatial streams

% V = 128x128 (Tx Tx) matrix of right singular vectors, represent the optimal "precoding" vectors for Tx
[U, D, V] = svd(H);

weightsTx = V(:,1:nLayers); % Calc transmitter weights (: = all rows, 1 upto nLayers columns)


weightsRx = U(:,1:nLayers)';

end

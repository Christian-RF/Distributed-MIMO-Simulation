function plotSitePatterns(tx, weightsPerBS, fc)
% PLOTSITEPATTERNS Plots Azimuth/Elevation cuts for each beam per BS.

numBS = numel(tx);
numLayers = size(weightsPerBS{1}, 2);

for b = 1:numBS
    figure('Name', sprintf('BS%d Beamforming Patterns', b), ...
           'Color', 'w', 'Position', [100+50*b, 100+50*b, 1200, 600]);
    sgtitle(sprintf('BS%d Tx Site: Beamforming Patterns (Beams 1-%d)', b, numLayers));
    
    for i = 1:numLayers
        % Set taper for current beam
        tx(b).Antenna.Taper = weightsPerBS{b}(:, i);
        
        % Top Row: Azimuth Cut
        subplot(2, numLayers, i);
        patternAzimuth(tx(b).Antenna, fc);
        title(sprintf('Beam %d Azimuth', i));
        
        % Bottom Row: Elevation Cut
        subplot(2, numLayers, i + numLayers);
        patternElevation(tx(b).Antenna, fc);
        title(sprintf('Beam %d Elevation', i));
    end
end
end
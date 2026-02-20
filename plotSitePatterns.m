function plotSitePatterns(txSuperCAntenna, txAcademicaAntenna, weightsSuperC, weightsAcademica, fc)
% PLOTSITEPATTERNS Plots Azimuth/Elevation cuts for 4 beams per site.
% Generates 2 Figures: one for SuperC, one for Academica.
% Each figure has a 2x4 layout (Top row: Azimuth, Bottom row: Elevation).

% Fig 1: SuperC
figure('Name', 'SuperC Beamforming Patterns', 'Color', 'w', 'Position', [100, 100, 1200, 600]);
sgtitle('SuperC Tx Site: Beamforming Patterns (Beams 1-4)');

for i = 1:4
    % Update Taper for the current beam index
    txSuperCAntenna.Taper = weightsSuperC(:, i);

    % Top Row: Azimuth Cut
    subplot(2, 4, i);
    patternAzimuth(txSuperCAntenna, fc);
    title(['Beam ' num2str(i) ' Azimuth']);

    % Bottom Row: Elevation Cut
    subplot(2, 4, i + 4);
    patternElevation(txSuperCAntenna, fc);
    title(['Beam ' num2str(i) ' Elevation']);
end

% Fig 2: Academica
figure('Name', 'Academica Beamforming Patterns', 'Color', 'w', 'Position', [150, 150, 1200, 600]);
sgtitle('Academica Tx Site: Beamforming Patterns (Beams 1-4)');

for i = 1:4
    % Update Taper for the current beam index
    txAcademicaAntenna.Taper = weightsAcademica(:, i);

    % Top Row: Azimuth Cut
    subplot(2, 4, i);
    patternAzimuth(txAcademicaAntenna, fc);
    title(['Beam ' num2str(i) ' Azimuth']);

    % Bottom Row: Elevation Cut
    subplot(2, 4, i + 4);
    patternElevation(txAcademicaAntenna, fc);
    title(['Beam ' num2str(i) ' Elevation']);
end
end
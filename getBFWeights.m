function W = getBFWeights(hEstimation, numLayers, scOffset, numRB)
%GETBEAMFORMINGWEIGHTS Erzeugt SVD-basierte Beamforminggewichte aus hest.
%
%   EINGABEN:
%   --------
%   hest      : 4D-Kanalarray der Form (subcarrier x symbol x rxAnt x txAnt),
%               erzeugt z.B. durch nrPerfectChannelEstimate.
%   nLayer    : Anzahl der gewünschten Datenströme (Layers).
%   scOffset  : Start-Subcarrier (Index, beginnt bei 0 oder 1, je nach Konvention).
%   noRBs     : Anzahl der RBs (Resource Blocks), über die der Kanal gemittelt wird.
%
%   AUSGABE:
%   -------
%   W         : Beamforminggewichte als (txAnt x nLayer)-Matrix, also die
%               Transmit-Precoder-Vektoren für nLayer Streams.
%
%   HINTERGRUND:
%   ------------
%   1) Extraktion/Aggregation des Kanals über best. Subcarrier-Bereich
%   2) Mittelung über alle (ausgewählten) Subcarrier und alle Symbole
%   3) SVD des resultierenden mittleren Kanals
%   4) Auswahl der ersten nLayer Singulärvektoren aus V als BF-Gewichte

    % --- 1. Festlegen, wie viele Subcarrier pro RB ---
    scPerRB = 12;   % Im 5G-NR-Standard meist 12 Subcarrier pro RB

    % --- 2. Bestimme Subcarrier-Indizes ---
    % Bei scOffset = 0 und noRBs = 1 würdest Du z.B. Subcarrier 1..12 nehmen.
    scIdx = (scOffset + 1) : (scOffset + numRB * scPerRB);

    % --- Sicherheits-Check: Grenzen nicht überschreiten ---
    maxSC = size(hEstimation,1);
    if any(scIdx > maxSC)
        error('Die gewünschten Subcarrier liegen außerhalb der Dimension von hest.');
    end

    % --- 3. Durchschnittliche Kanalantwort über diesen Frequenzbereich & alle Symbole ---
    % hest hat die Form (SC x Symbol x RxAnt x TxAnt).
    % => wir mitteln über SC und Symbol.
    Havg = mean(mean(hEstimation(scIdx, :, :, :), 1), 2);
    % Ergebnis hat nun die Form (1 x 1 x RxAnt x TxAnt).
    % => Dimension auf (RxAnt x TxAnt) reduzieren:
    Havg = squeeze(Havg);

    % --- 4. SVD auf mittlerem Kanal (RxAnt x TxAnt) ---
    % Havg = U * S * V'
    [U,S,V] = svd(Havg);

    % --- 5. Beamforming-Gewichte bestimmen: 
    %      Nehme die ersten nLayer Spalten von V (entspricht Tx-Seite).
    W = V(:, 1:numLayers);

    % OPTIONAL: Normierung der Spaltenvektoren (jeder Stream dieselbe Leistung)
    % for i = 1:nLayer
    %     W(:, i) = W(:, i) / norm(W(:, i));
    % end
end
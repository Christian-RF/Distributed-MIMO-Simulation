function [E_theta, E_phi] = getPolComponents(antenna, fc, az, el)
% Returns complex E-field components for each angle pair
% Output: [1 x numRays] each
    numRays = length(az);
    E_theta = zeros(1, numRays);
    E_phi   = zeros(1, numRays);
    
    for j = 1:numRays
        % step() on a polarized NRAntennaElement returns [2 x 1]:
        %   row 1 = H-pol (≈ phi component)
        %   row 2 = V-pol (≈ theta component)
        resp = step(antenna, fc, [az(j); el(j)]);
        
        % Map H/V to theta/phi (spherical coordinates)
        E_theta(j) = resp(2);   % V-pol → theta
        E_phi(j)   = resp(1);   % H-pol → phi
    end
end
function eff = efficiency_2D(angles, efficiency_metric)
%EFFICIENCY_2D Summary Calculates the SNR efficiency for a 3D radial acquisition
%compared to optimal sampling
%
%   
%
%   INPUTS:   angles    - list of angles of spokes (radians (zero to pi))
%             efficiency_metric - string. possible values:
%                                         'Winkelmann' (IEEE TRANSACTIONS ON MEDICAL IMAGING, 2007)
%                                         'electrostatic_potential'
%   OUTPUTS:  eff       - SNR efficiency of this sampling scheme compared 
%                         with an equal number of optimally spaced spokes
%
%   DEPENDENCIES: electrostatic_potential()
%
%
% Sophie Schauman, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(angles); % number of spokes
switch efficiency_metric
    case 'Winkelmann'
        spokes_relative_angles = angles-angles(1);
        spokes_unwrapped = mod(spokes_relative_angles,pi);

        q = diff(sort([spokes_unwrapped, pi])); %angle between adjacent spokes
        eff = sqrt((pi^2/N)/sum((0.5*(q+circshift(q,1))).^2));
        
    case 'electrostatic_potential'
        optimal_points = zeros(N, 2);
        tmp = linspace(0,0.5,N+1);
        optimal_points(:,2) = tmp(1:end-1);
        
        points = zeros(N, 2);
        spokes_relative_angles = angles-angles(1);
        spokes_unwrapped = mod(spokes_relative_angles,pi);
        points(:,2) = spokes_unwrapped/(2*pi);
        
        
        ref_pot = electrostatic_potential(optimal_points);
        eff = ref_pot/electrostatic_potential(points);
        
    otherwise
        error(['Efficiency metric `' efficiency_metric ...
            '` not implemented. --- Implemented metrics: '...
            '`Winkelmann`, `electrostatic_potential`'])
end
end









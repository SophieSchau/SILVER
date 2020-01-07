function eff = efficiency_3D(points, efficiency_metric)
%EFFICIENCY_3D Summary Calculates the SNR efficiency for a 3D radial acquisition
%compared to optimal sampling
%
%   
%
%   INPUTS:   points    - cartesian coordinates (x,y) of the distribution 
%                         of spokes onto a unit square.
%             efficiency_metric - string. possible values:
%                                         'voronoi_sphere'
%                                         'voronoi_plane'
%                                         'electrostatic_potential'
%   OUTPUTS:  eff       - SNR efficiency of this sampling scheme compared 
%                         with an equal number of optimally spaced spokes
%
%   DEPENDENCIES: voronoi_areas_on_sphere()
%                 voronoi_areas_unit_square()
%                 low_energy_configuration()
%                 electrostatic_potential()
%
%
% Sophie Schauman, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(points,1); % number of spokes
switch efficiency_metric
    case 'voronoi_sphere'
        points = uniquetol(points,1e-6, 'Byrows', true); % remove duplicates
        eff = 2*pi*sqrt(1./(N*sum((voronoi_areas_on_sphere(points).^2))));
        
    case 'voronoi_plane'
        points = uniquetol(points,1e-6, 'Byrows', true); % remove duplicates
        eff = sqrt(1./(N*sum((voronoi_areas_unit_square(points).^2))));
        
    case 'electrostatic_potential'
        ref_pot = low_energy_configuration(N);
        eff = ref_pot/electrostatic_potential(points);
        
    otherwise
        error(['Efficiency metric `' efficiency_metric ...
            '` not implemented. --- Implemented metrics: '...
            '`voronoi_sphere`, `voronoi_plane`, `electrostatic_potential`'])
end
end


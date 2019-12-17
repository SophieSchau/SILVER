function eff = efficiency_3D(points,N, efficiency_metric)
%EFFICIENCY_3D Summary Calculates the SNR efficiency for a 3D radial acquisition
%compared to optimal sampling
%
%   
%
%   INPUTS:   points    - cartesian coordinates (x,y) of the distribution 
%                         of spokes onto a unit square.
%             N         - number of spokes (number of points might be lower
%                         if there are duplicates)
%             efficiency_metric - string. possible values:
%                                         'voronoi_sphere'
%                                         'voronoi_plane'
%   OUTPUTS:  eff       - SNR efficiency of this sampling scheme compared 
%                         with an equal number of optimally spaced spokes
%
%   DEPENDENCIES: voronoi_areas_on_sphere()
%                 voronoi_areas_unit_square()
%
%
% Sophie Schauman, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch efficiency_metric
    case 'voronoi_sphere'
        eff = 2*pi*sqrt(1./(N*sum((voronoi_areas_on_sphere(points).^2))));
    case 'voronoi_plane'
        eff = sqrt(1./(N*sum((voronoi_areas_unit_square(points).^2))));
    otherwise
        error(['Efficiency metric `' efficiency_metric ...
            '` not implemented. --- Implemented metrics: '...
            '`voronoi_sphere`, `voronoi_plane`'])
end
end


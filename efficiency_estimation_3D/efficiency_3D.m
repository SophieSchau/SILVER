function eff = efficiency_3D(points,N)
%EFFICIENCY_3D Summary Calculates the SNR efficiency for a 3D radial acquisition
%compared to optimal sampling
%
%   
%
%   INPUTS:   points    - cartesian coordinates (x,y) of the distribution 
%                         of spokes onto a unit square.
%             N         - number of spokes (number of points might be lower
%                         if there are duplicates)
%   OUTPUTS:  eff       - SNR efficiency of this sampling scheme compared 
%                         with an equal number of optimally spaced spokes
%
%   DEPENDENCIES: voronoi_areas_on_sphere()
%
%
% Sophie Schauman, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    eff = 2*pi*sqrt(1./(N*sum((voronoi_areas_on_sphere(points).^2)))); 
end


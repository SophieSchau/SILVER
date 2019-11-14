function areas = voronoi_areas_on_sphere(points)
%VORONOI_AREAS_ON_SPHERE estimates the surface area occupide by each radial 
%spoke
%
%    Returns 2D area estimates for points based on voronoi parcellation on
%    sphere
% 
%    INPUTS:
%       points   -   Nx2 vector of points, scaled between [0,1]
%
%    OUTPUTS:
%       areas   -   Nx1 vector of densities for each point in 'points'
%
%    DEPENDENCIES:
%       'Voronoi Sphere' from https://uk.mathworks.com/matlabcentral/fileexchange/40989-voronoi-sphere
%
% Mark Chiew, Sophie Schauman 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of original points
N   = size(points,1);

% Transform points to points on unit sphere
q = zeros(N,3);
for i = 1:N
    q(i,3) = points(i,1);
    q(i,2) = sqrt(1-q(i,3).^2)*sin(2*pi*points(i,2));
    q(i,1) = sqrt(1-q(i,3).^2)*cos(2*pi*points(i,2));
end

% Compute voronoi cells
[~,~,~,areas] = voronoisphere(unique([q;-q],'rows')');
areas = areas(1:N);
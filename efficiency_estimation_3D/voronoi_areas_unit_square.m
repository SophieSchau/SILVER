function d = voronoi_areas_unit_square(p, display)
%VORONOI_AREAS_UNIT_SQUARE estimates the surface area occupide by each
%point, p, in the unit square
% Returns 2D density estimates for points based on voronoi parcellation
% 
% 	 INPUTS:
%       p   -   Nx2 vector of points, scaled between [0,1]
%       display - boolean. Show density distribution.
%
% 	 OUTPUTS
%       d   -   Nx1 vector of densities for each point in p
%
%
% Mark Chiew, Sophie Schauman 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    display = false;
end

% Number of original points
N   = size(p,1);
% Add points above, below, and diagonal to get cyclic boundary
pp  =   cat(1,p,...
            p+[-1, 0],...
            p+[+1, 0],...
            p+[ 0,-1],...
            p+[ 0,+1],...
            p+[-1,-1],...
            p+[-1,+1],...
            p+[+1,-1],...
            p+[+1,+1]);
% Compute voronoi cells
[v,c] = voronoin(pp);
% compte areas
d   =   zeros(N,1);
for i = 1:N
    d(i)    =   polyarea(v(c{i},1),v(c{i},2));
end
% Display density diagram if option set
if display
    figure();
    scatter(p(:,1),p(:,2),d*1E4,'filled');
end
end
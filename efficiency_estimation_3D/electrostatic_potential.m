function potential = electrostatic_potential(points)
%ELECTROSTATIC_POTENTIAL calculates the electrostatic potential of the
%points from a unit square transformed and duplicated onto each hemisphere
%of the unit sphere.
%
%   INPUTS:   points           - Nx2 array with x and y corrdinates of
%                                points on a unit square
%   OUTPUTS:  potential        - Total potential energy of the configuration
%
%
% Sophie Schauman, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert to points on hemisphere

w = zeros(size(points,1),3);
for i = 1:size(points,1)
    w(i,3) = points(i,1);
    w(i,2) = sqrt(1-w(i,3).^2)*sin(2*pi*points(i,2));
    w(i,1) = sqrt(1-w(i,3).^2)*cos(2*pi*points(i,2));
end
points = [w;-w]; %duplicate to get both ends of spoke

% Compute distance of every point from every other
% Use trick that ||x-y||^2 = ||x||^2 + ||y||^2 - 2x'y
d = abs(sqrt(2 - 2*points*points'));

% Sum over upper right triangle
% This accounts for every pair-wise interaction
potential = 0;
for i = 1:size(points,1)-1    
    potential = potential + sum(1./d(i,i+1:end));
end
if potential == Inf
    potential = realmax;
end
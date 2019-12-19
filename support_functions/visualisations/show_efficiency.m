function [] = show_efficiency(ratio, window_sizes, efficiency_metric) 
%SHOW_EFFICIENCY plots efficiency as a function of window size for radial
%spokes with set increments in 2D or 3D
% 
%    INPUTS:
%       ratio   -   a single number (2D imaging) or an array of two numbers
%                   (3D imaging). For 3D, the first number is the ratio
%                   that determines the azimuthal angle of the spoke, and 
%                   the second determines z-coordinate of spoke tip.
%       window_sizes  - list of window sizes to calculate efficiency for
%                   and plot.
%       efficiency_metric - for 3D imaging this needs to be defined as
%                   there is no universal metric. Can be left empty for 2D.
%                   'voronoi_sphere' is default.
%    DEPENDENCIES:
%       efficiency_range()
%
% Sophie Schauman 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    efficiency_metric = 'voronoi_sphere';
end

eff_a = efficiency_range(ratio,window_sizes,efficiency_metric) ;
    plot(window_sizes, eff_a, 'Linewidth',3)
end
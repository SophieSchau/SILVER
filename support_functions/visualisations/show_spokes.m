function [] = show_spokes(ratio, spoke_numbers) 
%SHOW_SPOKES plots spokes generated by a set increment method (e.g. SILVER 
%or Golden means).
% 
%    INPUTS:
%       ratio   -   a single number (2D imaging) or an array of two numbers
%                   (3D imaging). For 3D, the first number is the ratio
%                   that determines the z-coordinate of spoke tip, and 
%                   the second determines azimuthal angle of the spoke.
%       spoke_numbers  - the indecies of spokes to plot
%
% Sophie Schauman 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax = gca;
linecolor = ax.ColorOrder(ax.ColorOrderIndex,:);
colorindex = ax.ColorOrderIndex;
if length(ratio) == 1
    angles = mod(spoke_numbers*ratio,1)*pi;
    x = cos(angles);
    y = sin(angles);
    h = plot([x;-x], [y;-y], 'Linewidth', 2, 'color', linecolor);
elseif length(ratio) == 2
    error('3D not implemented yet')
end
ax.ColorOrderIndex = mod(colorindex,7)+1;
set(h(2:end), 'Handlevisibility', 'off')
end
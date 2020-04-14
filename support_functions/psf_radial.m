function [psf_map] = psf_radial(ratio, N, matrix_size, sens)
%PSF_RADIAL generates a point spread function map for radially sampled MRI
%   
%   INPUTS: ratio        -   angular step size is ratio*180 degrees.
%           N            -   number of spokes to generate
%           matrix_size  -   size of psf_map. Defaults to [N*pi/2 x N*pi/2]
%           sens         -   coil sensityivity map. Defaults to single coil
%
%   OUTPUT: psf_map -   2D map of point spread function. 
%
% Sophie Schauman 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 3
    matrix_size = [round(N*pi/2), round(N*pi/2)];
end
if nargin < 4
    sens = ones(matrix_size);
end


if length(matrix_size) == 1
    matrix_size = [matrix_size,matrix_size];
end
point_map = zeros(matrix_size);
point_map(round(matrix_size(1)./2),round(matrix_size(2)./2)) = 1;

angles = (0:N-1)*ratio*pi;

kspace = gen_radial_traj(angles, max(matrix_size)*2, []);
kspace= reshape(kspace,[],1,2);

E = xfm_NUFFT([size(point_map) 1 1],[],[],kspace,'wi', 1);

psf_map = E'.*(E.*point_map);
    


end


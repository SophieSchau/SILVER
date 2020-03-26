function [psf_map] = psf_radial(ratio,N, matrix_size)
%PSF_RADIAL generates a point spread function map for radially sampled MRI
%   
%   INPUTS: ratio        -   angular step size is ratio*180 degrees.
%           N            -   number of spokes to generate
%           matrix_size  -   size of psf_map. Defaults to [N*pi/2 x N*pi/2]
%
%   OUTPUT: psf_map -   2D map of point spread function. 
%
% Sophie Schauman 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 3
    matrix_size = [round(N*pi/2), round(N*pi/2)];
end
point_map = zeros(matrix_size);
point_map(round(matrix_size./2)) = 1;

angles = mod((0:N-1)*ratio,1)*pi;

%kspace = gen_radial_traj(angles, N, []);
%kspace= reshape(kspace,[],1,2);

E = xfm_NUFFT([round(N*pi/2) round(N*pi/2) 1 1],[],[],kspace);

psf_map = E'.*(E.*point_map);
    






end


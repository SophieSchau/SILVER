function [psf] = psf_radial(ratio, N, matrix_size)
%PSF_RADIAL generates a point spread function map for radially sampled MRI
%(single coil)
%   
%   INPUTS: ratio        -   angular step size is ratio*180 degrees.
%           N            -   number of spokes to generate
%           (matrix_size)-   Optional size of matrix to reconstruct.
%                            Defaults to round(sqrt(2)*N/pi).
%
%   OUTPUT: psf -   2D matrix. Operator that when applied to a delta
%                   function characterizes the effect of the chosen
%                   trajectory has on a point source.
%
% Sophie Schauman 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 3
    matrix_size = round(sqrt(2)*N/pi);
end

angles = (0:N-1)*ratio*pi;

k = sqrt(2)*gen_radial_traj(angles, max(matrix_size)*2, []);

[x, y]  = meshgrid(-floor(matrix_size/2):ceil(matrix_size/2)-1, -floor(matrix_size/2):ceil(matrix_size/2)-1);
E       = exp(-1j*(k(:,2)*x(:)' + k(:,1)*y(:)'));
psf     = E'*E;

end


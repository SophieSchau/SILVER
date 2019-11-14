function eff = efficiency_2D(spokes)
%EFFICIENCY_2D Calculates the SNR efficiency for a 2D radial acquisition
%compared to uniform sampling
%
%   Based on Winkelmann et al. (IEEE TRANSACTIONS ON MEDICAL IMAGING, 2007)
%
%   INPUTS:   spokes - list of angles of spokes (radians)
%   OUTPUTS:  eff - SNR efficiency of this sampling scheme compared with an
%                   equal number of uniformly spaced spokes
%
%
% Sophie Schauman, Mark Chiew, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    spokes_relative_angles = spokes-spokes(1);
    spokes_unwrapped = mod(spokes_relative_angles,pi);

    q = diff(sort([spokes_unwrapped, pi])); %angle between adjacent spokes
    eff = sqrt((pi^2/length(spokes))/sum((0.5*(q+circshift(q,1))).^2));
end
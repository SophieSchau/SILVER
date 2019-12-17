
function effs = efficiency_range(ratio,window_sizes) 
%EFICIENCY_RANGE get a list of SNR efficiencies for each number of spokes
%in the range [M, N]
%
%   INPUTS:   ratio - a single number for 2D acquisition (equivalent to
%             golden ratio), a 2x1 vector for 3D acquisition (equivalen to
%             the 2D Golden means (Chan, 2009))
%             window_sizes - list of number of spokes in temporal window
%   OUTPUTS:  effs - list of efficiencies for window size, M to N.
%
%   DEPENDENCIES: - efficiency_2D()
%                 - surface_points()
%
% Sophie Schauman, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    effs = zeros(length(window_sizes),1);
    
    if length(ratio) == 1 % 2D acquisition
    
        for ii = 1:length(effs)
            effs(ii) = efficiency_2D([0:window_sizes(ii)-1]*ratio*pi);        
        end   

    
    elseif length(ratio) == 2 % 3D acquisition
        for ii = 1:length(effs)
            
            x = mod([0:window_sizes(ii)-1]*ratio(1),1)';
            y = mod([0:window_sizes(ii)-1]*ratio(2),1)';
            points = [x,y]; % tip positions on unit square
            points = uniquetol(points,1e-6, 'Byrows', true);
            effs(ii) = efficiency_3D(points, window_sizes(ii)); 
        end   
    end
end
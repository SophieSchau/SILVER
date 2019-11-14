
function effs = efficiency_range(ratio,M,N) 
%EFICIENCY_RANGE get a list of SNR efficiencies for each number of spokes
%in the range [M, N]
%
%   INPUTS:   ratio - a single number for 2D acquisition (equivalent to
%             golden ratio), a 2x1 vector for 3D acquisition (equivalen to
%             the 2D Golden means (Chan, 2009))
%             M - minimum number of spokes in temporal window
%             N - maximum number of spokes in temporal window
%   OUTPUTS:  effs - list of efficiencies for window size, M to N.
%
%   DEPENDENCIES: - efficiency_2D()
%
% Sophie Schauman, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    effs = zeros(N-M+1,1);
    
    if length(ratio) == 1 % 2D acquisition
    
        for i = 1:length(effs)
            effs(i) = efficiency_2D((0:M-2+i)*ratio*pi);        
        end   

    
    elseif length(ratio) == 2 % 3D acquisition
        error('3D acquisition not implemented yet')
    end
end
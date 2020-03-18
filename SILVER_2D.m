function [ratio] = SILVER_2D(window_sizes, efficiency_metric, savefilename)
%SILVER_2D Calculate the optimal increment for a range of temporal windows
%   For Golden angle radial sampling the optimal angular increment between 
%   spokes is pi times the Golden ratio (approximately 0.6180). The SILVER
%   method optimises the increment for a window size [M,N]. 
%   
%   The optimisation is performed by a minimized of the maximum SNR 
%   inefficiency (as the reciprocal of the SNR efficiency ?). 
%
%   INPUTS:   window_sizes - list of number of spokes used to reconstruct
%             efficiency_metric - string. possible values:
%                                         'Winkelmann' (default)
%                                         'electrostatic_potential'
%             savefilename - name and path of file to save result to.
%   OUTPUTS:  ratio - the SILVER equivalent of the golden ratio. The
%                     optimal angle increment is pi * ratio.
%
%   DEPENDENCIES: - efficiency_range()
%                   gr2D()
%
% Sophie Schauman, Mark Chiew, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    opts = optimoptions('fmincon','Display','none','Algorithm','interior-point');    
    ratios = zeros(100,1);
    fval = zeros(100,1);
    
    
    % Loop starting estimates to overcome non-convexity 
    % This should be paralellised for efficiency
    rng(1)
    tic
    for i = 1:100
        if i == 1
            [ratios(i), fval(i)] = fmincon(@(x)1./min(efficiency_range(x,window_sizes, efficiency_metric)),gr2D, [], [], [], [], 0, 1, [], opts);        
        else
            [ratios(i), fval(i)] = fmincon(@(x)1./min(efficiency_range(x,window_sizes, efficiency_metric)),rand(1), [], [], [], [], 0, 1, [], opts);        
        end 
            
    end
    t = toc;
    [~,i] = min(fval);
    ratio = ratios(i);
    eff_Golden = efficiency_range(gr2D,window_sizes,efficiency_metric);
    eff_SILVER  = efficiency_range(ratio,window_sizes,efficiency_metric);
    
    fprintf(1,'Golden Angle: %g, Min SNR-Efficiency: %g\n',180*gr2D,min(eff_Golden));
    fprintf(1,'SILVER Angle: %g, Min SNR-Efficiency: %g\n',rad2deg(ratio*pi),min(eff_SILVER));
    
    if exist('savefilename', 'var')
        save(savefilename, 'eff_SILVER', 'fval', 'ratio', 'ratios', 't', 'window_sizes')
    end

end



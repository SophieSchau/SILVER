function [ratio] = SILVER_2D(window_sizes, efficiency_metric, savefilename, max_spokes, mat_size)
%SILVER_2D Calculate the optimal increment for a range of temporal windows
%   For Golden angle radial sampling the optimal angular increment between 
%   spokes is pi times the Golden ratio (approximately 0.6180). The SILVER
%   method optimises the increment for a set of window sizes, S. 
%   
%   The optimisation is performed by a minimization of the maximum SNR 
%   inefficiency (as the reciprocal of the SNR efficiency ?). 
%
%   This function also includes a check that the Nyquist criterion is met
%   when combining at most max_spokes number of spokes. This is important
%   for time-avergared coil combination estimation.
%
%   INPUTS:   window_sizes - list of number of spokes used to reconstruct
%             efficiency_metric - string. possible values:
%                                         'Winkelmann' (default)
%                                         'electrostatic_potential'
%             savefilename - name and path of file to save result to.
%             max_spokes - maximum number of spokes that the nyquist limit
%                          has to be reached within.
%             mat_size - maximum dimension of matrix size to be 
%                        reconstructed. This is for calculation of how many
%                        uniformly distributed spokes are needed for the
%                        nyquist limit (pi/2*mat_size).
%   OUTPUTS:  ratio - the SILVER equivalent of the golden ratio. The
%                     optimal angle increment is pi * ratio.
%
%   DEPENDENCIES: - efficiency_range()
%
% Sophie Schauman, Mark Chiew, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 5
        mat_size = [];
    end
    if nargin < 4
        max_spokes = [];
    end
    if nargin < 3
        savefilename = [];
    end
    
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
    
    % check for nyquist limit
    if ~isempty(max_spokes)
        % find out maximum gap
        nyquist_nSpokes = pi/2*mat_size; % for uniform radial sampling (Bernstein)
        nyquist_gapSize = (1/nyquist_nSpokes*pi)*1.01; % radians - allow for 1% deviation

        angles = [0:max_spokes-1]*ratio*pi;
        spokes_relative_angles = angles-angles(1);
        spokes_unwrapped = mod(spokes_relative_angles,pi);

        q = diff(sort([spokes_unwrapped, pi]));
        max_gapSize = max(q);
        
        while max_gapSize > nyquist_gapSize
            [~,i] = min(fval);
            fval(i) = inf; % make this result ineligible
            [~,i] = min(fval); % choose next best option
            
            ratio = ratios(i);
            
            angles = [0:max_spokes-1]*ratio*pi;
            spokes_relative_angles = angles-angles(1);
            spokes_unwrapped = mod(spokes_relative_angles,pi);

            q = diff(sort([spokes_unwrapped, pi]));
            max_gapSize = max(q);
            
            if min(fval) == inf
                ratio = [];
                warning('No step size that could reach the nyquist limit within the maximum number of spokes could be found')
                break
            end
        end

    
%     eff_Golden = efficiency_range(gr2D,window_sizes,efficiency_metric);
%     eff_SILVER  = efficiency_range(ratio,window_sizes,efficiency_metric);
%     
%     fprintf(1,'Golden Angle: %g, Min SNR-Efficiency: %g\n',180*gr2D,min(eff_Golden));
%     fprintf(1,'SILVER Angle: %g, Min SNR-Efficiency: %g\n',rad2deg(ratio*pi),min(eff_SILVER));
    
    if ~isempty(savefilename)
        save(savefilename, 'eff_SILVER', 'fval', 'ratio', 'ratios', 't', 'window_sizes')
    end

end



function [ratio, spokes_for_nyquist] = SILVER_2D_nyquist(window_sizes, efficiency_metric, max_spokes, mat_size, savefilename)
%SILVER_2D_nyquist Calculate the optimal increment for a range of temporal windows
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
%             max_spokes - maximum number of spokes that the nyquist limit
%                          has to be reached within.
%             mat_size - maximum dimension of matrix size to be 
%                        reconstructed. This is for calculation of how many
%                        uniformly distributed spokes are needed for the
%                        nyquist limit (pi/2*mat_size).
%             savefilename - name and path of file to save result to.
%   OUTPUTS:  ratio - the SILVER equivalent of the golden ratio. The
%                     optimal angle increment is pi * ratio.
%             spokes_for_nyquist - how many SILVER ordered spokes are
%                                  needed to reach the nyquist limit. This 
%                                  number must be smaller or equal to 
%                                  max_spokes.
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
for i = 1:100 % 100 different starting values. One of them is gr2D
    if i == 1
        [ratios(i), fval(i)] = fmincon(@(x)1./min(efficiency_range(x,window_sizes, efficiency_metric)),gr2D, [], [], [], [], 0, 1, [], opts);        
    else
        [ratios(i), fval(i)] = fmincon(@(x)1./min(efficiency_range(x,window_sizes, efficiency_metric)),rand(1), [], [], [], [], 0, 1, [], opts);        
    end 

end
t = toc;


acceptable_result = false;
nyquist_limit = pi/2*mat_size; % for uniform radial sampling (Bernstein)
gap_nyquist = (1/nyquist_limit*pi)*1.01; % radians - allow for 1% deviation


for N = 1:100 % check at most 100 different results of the optimisation
    if acceptable_result
        break
    else
        [~,i] = min(fval);
        ratio = ratios(i);
        eff_Golden = efficiency_range(gr2D,window_sizes,efficiency_metric);
        eff_SILVER  = efficiency_range(ratio,window_sizes,efficiency_metric);

        % How many spokes are needed such that the nyquist criterion is
        % reached?

        for Nspokes = 2:max_spokes
            angles = [0:Nspokes-1]*ratio*pi;
            spokes_relative_angles = angles-angles(1);
            spokes_unwrapped = mod(spokes_relative_angles,pi);

            q = diff(sort([spokes_unwrapped, pi]));
            max_gap = max(q);
            
            if max_gap <= gap_nyquist
                acceptable_result = true;
                spokes_for_nyquist = Nspokes;
                fprintf(1,'Golden Angle: %g, Min SNR-Efficiency: %g\n',180*gr2D,min(eff_Golden));
                fprintf(1,'SILVER Angle: %g, Min SNR-Efficiency: %g\n',rad2deg(ratio*pi),min(eff_SILVER));
    
                break
            end
        end
        
        if max_gap > gap_nyquist
            acceptable_result = false;
            fval(i) = inf; %choose next best result
        end
    end

    
    
if exist('savefilename', 'var')
    save(savefilename, 'eff_SILVER', 'fval', 'ratio', 'ratios', 't', 'window_sizes')
end

end

if ~acceptable_result
    ratio = 0;
    disp('No result which reached the nyquist criterion could be found')
end


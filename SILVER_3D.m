function [ratio, eff_gr,eff_a] = SILVER_3D(window_sizes, efficiency_metric)
%SILVER_3D Calculate the optimal increment for a range of temporal windows
%   For 2D Golden means radial sampling the azimuthal angle, and and 
%   z-increment between the tip position of subsequently acquired spokes is
%   determined by the 2D Golden means (Chan, 2009):
%
%   phi_n = phi_n-1 + 2*pi*mod(0.6823,1)
%   z_n = mod(z_n-1 +0.4656,1)
%
%   The SILVER method optimises the increment for a window size [M,N] 
%   instead of [1 to infinity]. 
%   
%   The optimisation is performed by a minimized of the maximum SNR 
%   inefficiency (as the reciprocal of the SNR efficiency ?). 
%
%   INPUTS:   window_sizes - list of number of spokes used to reconstruct
%             efficiency_metric - string. possible values:
%                                         'voronoi_sphere' (default)
%                                         'voronoi_plane'
%   OUTPUTS:  ratio - the SILVER equivalent of the 2D golden means. 
%                     (a 2x1 vector) 
%
%   DEPENDENCIES: - efficiency_range()
%                   gr3D()
%
% Sophie Schauman, Mark Chiew, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% a is the optimal angle increment for an arbitrary number of spokes
% between N-M for a 3D acquisition

    if nargin < 2
        efficiency_metric = 'voronoi_sphere';
    end
    disp('3D acquisition')
    disp(['Efficiency defined by: ' efficiency_metric]);
    opts = optimoptions('fmincon','Display','none','Algorithm','interior-point');    
    ratios = zeros(2,25);
    fval = zeros(1,25);
    % Loop starting estimates to overcome non-convexity
    % should be parallelised for efficiency
    n = 1;
    for ii = 1:2:10
        for jj = 2:2:10
            [ratios(:,n),fval(n)] = fmincon(@(x)1./min(efficiency_range(x,window_sizes,efficiency_metric)),[(ii/10); (jj/10)], [], [], [], [], [0,0], [1,1], [], opts); 
            disp(['iteration: ' num2str(n) ' of 25'])
            n = n+1;
        end
    end
    
    [~,i] = min(fval);
    ratio = ratios(:,i);
    eff_gr = efficiency_range(gr3D,window_sizes,efficiency_metric);
    eff_a  = efficiency_range(ratio,window_sizes,efficiency_metric);
    
    gr = gr3D;
    
    angle_gr = gr(2)*2*pi;
    z_ratio_gr = gr(1);
    
    
    angle = ratio(2)*2*pi;
    z_ratio = ratio(1);
    
    fprintf(1,'Golden - Azimuthal Angle: %g, Z increment: %g, Min SNR-Efficiency: %g\n',rad2deg(angle_gr), z_ratio_gr,min(eff_gr));
    fprintf(1,'Silver - Azimuthal Angle: %g, Z increment: %g, Min SNR-Efficiency: %g\n',rad2deg(angle), z_ratio,min(eff_a));
    
    

end



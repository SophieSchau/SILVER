function [prcnt_increase, min_eff_SILVER, min_eff_GR, SILVER_ratios,b] = SILVER_vs_GR_discrete_sets(S, efficiency_metric, savename)
%SILVER_VS_GR_DISCRETE_SETS Does SILVER perform better than the golden 
%   ratio when the set of window sizes to optimize for are discrete 
%   temporal resolutions?
%                                                 
%   S can be a set or multiple sets of window sizes (cell array of
%   arrays).
%   
%   savename is used to save the result OR read a pre-calculated result.
                                             
%   Sophie Schauman 2020                       
%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    savename = [];
end

% Do SILVER optimization for the sets of ranges
n = 0;
for s = S
   n = n+1;
   savename_ratio = ['experiments/precalculated_' efficiency_metric '/silver_' strrep(num2str(s{:}),' ', '_') '.mat'];
   if ~exist(savename_ratio, 'file')
        SILVER_ratios(n) = SILVER_2D(s{:},efficiency_metric, savename_ratio) ;
   end
   disp([num2str(s{:}) ' done'])
end


if ~exist([savename '.mat'], 'file')
    n = 0;
    for s = S
        n = n+1;
        load(['experiments/precalculated_' efficiency_metric '/silver_' strrep(num2str(s{:}),' ', '_') '.mat'],'eff_SILVER','ratio')
            min_eff_SILVER(n) = (min(eff_SILVER));
            min_eff_GR(n) = min(efficiency_range(gr2D,s{:},'electrostatic_potential'));
            prcnt_increase(n) = ((min_eff_SILVER(n)/min_eff_GR(n))-1)*100;
            SILVER_ratios(n) = ratio;
    end

    if ~isempty(savename)
        save(savename, 'prcnt_increase','SILVER_ratios','min_eff_GR', 'min_eff_SILVER','S')
    end
else
    warning('using pre-calculated results')
    load(savename);
end
 




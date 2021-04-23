function [prcnt_increase, min_eff_SILVER, min_eff_GR, SILVER_ratios] = SILVER_vs_GR_ranges(M, R, efficiency_metric, savename)
%SILVER_VS_GR_RANGES For what ranges of window sizes does SILVER perform 
%better than the golden ratio method?
%
%   Starting from a minimum window size, M, how does the performance of
%   SILVER compare to GR (based on some metric) when the set included in 
%   the optimisation grows from S = {M and M+1} to all intermediate window 
%   sizes between M and M+R?
%
%   In this function M can be a single starting minimum or an array of many
%   starting minima. R is an array containing all set lengths.
%   
%   savename is used to save the result OR read a pre-calculated result.
                                             
%   Sophie Schauman 2020                       
%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    savename = [];
end

% Do SILVER optimization for the sets of ranges
for m = M
    for r = R
       savename_ratio = ['experiments/precalculated_' efficiency_metric '/silver_' num2str(m) 'to' num2str(m+r) '.mat'];
       if ~exist(savename_ratio, 'file')
            SILVER_2D(m:m+r,efficiency_metric, savename_ratio) ;
       end
       disp([num2str(m) 'to' num2str(m+r) ' done'])
    end
end


if ~exist([savename '.mat'], 'file')
    for w = 1:length(M)
        m = M(w);
        for r = R
            load(['experiments/precalculated_' efficiency_metric '/silver_' num2str(m) 'to' num2str(m+r) '.mat'],'eff_SILVER','ratio')
            min_eff_SILVER(w,r) = (min(eff_SILVER));
            min_eff_GR(w,r) = min(efficiency_range(gr2D,m:m+r,'electrostatic_potential'));
            prcnt_increase(w,r) = ((min_eff_SILVER(w,r)/min_eff_GR(w,r))-1)*100;
            SILVER_ratio(w,r) = ratio;
        end
    end
    if ~isempty(savename)
        save(savename, 'prcnt_increase','SILVER_ratio','min_eff_GR', 'min_eff_SILVER','M', 'R')
    end
else
    load(savename, 'prcnt_increase','SILVER_ratio','min_eff_GR', 'min_eff_SILVER','M', 'R')
    warning('using pre-calculated results')
end
    




figure
if length(R) > 10
    tickloc = [R(1) R(10:10:end)];
else
    tickloc = R;
end
xticks(tickloc)

box on
ii = 0;
for r = tickloc
    ii = ii+1;
    labels{ii} = ['S = \{M, ..., M+' num2str(r) '\}'];
end
set(gca,'xticklabel',labels)
xtickangle(-90)
hold on

for i = 1:length(M)
    c_map = i./length(M)*[0.7 0.7 0.7; 1 0.5 0];
    plot(R,min_eff_SILVER(i,R), 'o-','Linewidth', 3, 'markersize', 5, 'Color',c_map(1,:), 'DisplayName', ['SILVER, M = ' num2str(M(i))])
    hold on
    plot(R,min_eff_GR(i,R), '--','Linewidth', 2, 'markersize', 5, 'Color',[c_map(2,:)], 'DisplayName', ['Golden ratio, M = ' num2str(M(i))])
end



set(gca,'FontSize',16)
set(gca, 'LineWidth', 2)
grid on
set(gcf,'Position',[124 357 876 441])
axis([min(R)-1,max(R)+1,ylim])
legend
ylabel('Minimum efficiency, \eta')

if ~isempty(savename)
    savefig([savename '.fig'])
    saveas(gcf,[savename '.svg'])
end


function [] = SILVER_vs_GR_optimisation_example_fig(S, efficiency_metric, xrange, ylims, savename)
%SILVER_VS_GR_OPTIMISATION_EXAMPLE_FIG How does efficiency vary within a set of
%   window sizes ofor SILVER and for GR?
%
%   S can be a set or multiple sets of window sizes (cell array of
%   arrays).
%   xrange is the range of window sizes to visualise. ylims is the y-axis
%   limits.
%   
%   savename is used to save the result OR read a pre-calculated result.
                                             
%   Sophie Schauman 2020                       
%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    savename = [];
end

% Do SILVER optimization for the sets of ranges
for s = S
   savename_ratio = ['experiments/precalculated_' efficiency_metric '/silver_' strrep(num2str(s{:}),' ', '_') '.mat'];
   if ~exist(savename_ratio, 'file')
        SILVER_2D(s{:},efficiency_metric, savename_ratio) ;
   end
   disp([num2str(s{:}) ' done'])
end


eff_GR = efficiency_range(gr2D,xrange,efficiency_metric);
hold on

n = 0;
lg={};
for s = S
    n = n+1;
    c_map = n./length(S)*[0.7 0.7 0.7];
    for m = 1:length(s{:})
        rectangle('Position',[s{:}(m)-0.5,0.9,1,0.1], 'EdgeColor', 'None', 'CreateFcn', @(l, e) set(l, 'FaceColor', [c_map, 0.3]))
    end
    
    
    load(['experiments/precalculated_' efficiency_metric '/silver_' strrep(num2str(s{:}),' ', '_') '.mat'], 'ratio','eff_SILVER')
    eff_SILVER_all = efficiency_range(ratio,xrange,'electrostatic_potential');
    h = plot(xrange,eff_SILVER_all,'x-','Linewidth', 3, 'markersize', 5, 'color',c_map );
    lg{n} = ['SILVER: S = \{' num2str(s{:}) '\}, \alpha = ' num2str(ratio)];
    plot(s{:},eff_SILVER,'o','markersize', 10, 'markerfacecolor',c_map,'markeredgecolor',c_map,'HandleVisibility','off');

end
plot(xrange,eff_GR,'x-','Linewidth', 3, 'markersize', 5, 'color', [1 0.5 0])
lg{end+1} = ['Golden ratio: \alpha = ' num2str(gr2D)];

box on
grid minor
set(gcf,'Position',[124 359 876 439])
legend(lg{:}, 'location', 'southeast')
xticks(xrange);
xlabel('Window size')
ylabel('Efficiency, \eta')
axis([xlim,ylims])

if ~isempty(savename)
    savefig([savename '.fig'])
    saveas(gcf,[savename '.svg'])
end


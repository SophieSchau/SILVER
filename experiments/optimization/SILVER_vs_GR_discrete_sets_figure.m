function [idx_sc] = SILVER_vs_GR_discrete_sets_figure(S, data_filename, savename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


load(data_filename);


[~, idx_sc] = sort(prcnt_increase,'descend');
prcnt_increase = prcnt_increase(idx_sc);

c_map = [0.5 0.5 0.5; 1 0.5 0];
box on

yyaxis right

for n = 1:length(S)
    labels{n} = ['S = \{' num2str(S{n}) '\}'];
end
xtick([1:length(S)])
set(gca,'xticklabel',labels)
xtickangle(-90)
hold on

b = bar([min_eff_SILVER(idx_sc);min_eff_GR(idx_sc)]','barwidth', 1);
b(1).FaceColor = c_map(1,:);
b(2).FaceColor = c_map(2,:);
b(1).FaceAlpha = 0.5;
b(2).FaceAlpha = 0.5;
b(1).LineWidth = 0.5;
b(2).LineWidth = 0.5;

set(gca,'xticklabel',labels(idx_sc))
xtickangle(-90)
legend('SILVER', 'Golden ratio','location', 'northeast')

axis([xlim, 0.9,1.01])
grid minor
set(gcf,'Position',[124 357 876 441])

ylabel('Minimum efficiency, \eta')
set(gca, 'ycolor', 'k')




yyaxis('left')
b2 = bar([prcnt_increase'], 'EdgeColor', 'k', 'LineWidth', 1);
ylabel('% improvement in \eta')
set(gca, 'ycolor', 'k')
set(gca, 'SortMethod', 'depth')
ylim([1,3.5].*ylim);

big_increase = prcnt_increase>1;



for n = 1:length(S)
    if prod(S{n}(2:end)-S{n}(1:end-1))~=1
        labels{n}= ['S = \{' num2str(S{n}) '\}'];
    else
        labels{n} = ['S = \{' num2str(S{n}(1)) ' to ' num2str(S{n}(end)) '\}'];
    end
end


labels = labels(idx_sc);
xticklabels(labels)
yl = ylim;
for n = 1:length(big_increase)

        text(n,b2.YData(n)+0.35,[num2str(prcnt_increase(n),2) '%'], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold')


end

line(find(big_increase, 1,'last')+0.5.*[1,1], ylim, 'linewidth', 1, 'color', 'k', 'linestyle', '--', 'handlevisibility', 'off')

ll = legend;
ll.String{1} = 'Difference (left y-axis)';
ll.String{2} = 'SILVER (right y-axis)';
ll.String{3} = 'GR (right y-axis)';

set(ll, 'location', 'northeast')
set(gcf, 'position', [16 422 967 383]);


if ~isempty(savename)
    savefig([savename '.fig'])
    saveas(gcf,[savename '.svg'])
end

end


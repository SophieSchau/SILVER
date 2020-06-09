%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   For what ranges of window sizes does SILVER perform better than the 
%   golden ratio method?
%                                              
%   Sophie Schauman 2020                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% 1. Choose minimum window size, M, and length of ranges, R, to try

M = [4,16,32];
R = 1:100;

%% 2. Do SILVER optimization for those ranges, electrostatic potential cost
for m = M
    for r = R
       savename = ['examples/precalculated/silver_' num2str(m) 'to' num2str(m+r) '.mat'];
       if ~exist(savename, 'file')
            SILVER_2D(m:m+r,'electrostatic_potential', savename) ;
       end
       disp([num2str(m) 'to' num2str(m+r) ' done'])
    end
end

%% 3. Calculate minimum efficiency in ranges for SILVER and golden ratio

savename= 'examples/example1_ranges/example1_ranges_result.mat';

if ~exist(savename, 'file')
    for w = 1:length(M)
        m = M(w);
        for r = R
            load(['examples/precalculated/silver_' num2str(m) 'to' num2str(m+r) '.mat'],'eff_SILVER','ratio')
            min_eff_SILVER(w,r) = (min(eff_SILVER));
            min_eff_GR(w,r) = min(efficiency_range(gr2D,m:m+r,'electrostatic_potential'));
            prcnt_increase(w,r) = ((min_eff_SILVER(w,r)/min_eff_GR(w,r))-1)*100;
            SILVER_ratio(w,r) = ratio;

        end
    end
    save(savename, 'prcnt_increase','SILVER_ratio','min_eff_GR', 'min_eff_SILVER','M', 'L')
else
    warning('using pre-calculated results')
    load(savename);
end

%% 4. Vizualise result
figure(1)
xticks([1,10:10:100])
box on
ii = 0;
for n = [1,10:10:100]
    ii = ii+1;
    labels{ii} = ['S = \{M, ..., M+' num2str(R(n)) '\}'];
end
set(gca,'xticklabel',labels)
xtickangle(-90)
hold on

c_map = [0.7 0.7 0.7; 0.5 0.5 0.5; 0.3 0.3 0.3; 1 0.5 0;0.7*([1 0.5 0]); 0.5*([1 0.5 0])];
plot(R,min_eff_SILVER(1,:), 'o-','Linewidth', 3, 'markersize', 5, 'Color',c_map(1,:))
hold on
plot(R,min_eff_GR(1,:), '--','Linewidth', 2, 'markersize', 5, 'Color',[c_map(4,:)])
plot(R,min_eff_SILVER(2,:),'o-','Linewidth', 3,'markersize', 5, 'Color',c_map(2,:))
plot(R,min_eff_GR(2,:), '--','Linewidth', 2, 'markersize', 5, 'Color',[c_map(5,:)])
plot(R,min_eff_SILVER(3,:),'o-', 'Linewidth', 3,'markersize', 5, 'Color',c_map(3,:))
plot(R,min_eff_GR(3,:), '--','Linewidth', 2, 'markersize', 5, 'Color',[c_map(6,:)])


set(gca,'FontSize',16)
set(gca, 'LineWidth', 2)
grid on
set(gcf,'Position',[124 357 876 441])
axis([0,max(R)+1,ylim])
legend('SILVER, M = 4', 'Golden ratio, M = 4', 'SILVER, M = 16', 'Golden ratio, M = 16', 'SILVER, M = 32', 'Golden ratio, M = 32','location', 'northeast')
ylabel('Minimum efficiency, \eta')


savefig('examples/example1_ranges/example1_ranges_result.fig')
saveas(gcf,'examples/example1_ranges/example1_ranges_result.tiff')


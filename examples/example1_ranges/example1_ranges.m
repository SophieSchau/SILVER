%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   For what ranges of window sizes does SILVER perform better than the 
%   golden ratio method?
%                                              
%   Sophie Schauman 2020                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% 1. Choose minimum window size, M, and length of ranges, L, to try

M = [4,16,32];
L = 1:100;

%% 2. Do SILVER optimization for those ranges
for m = M
    for l = L
       savename = ['examples/precalculated/silver_' num2str(m) 'to' num2str(m+l) '.mat'];
       if ~exist(savename, 'file')
            SILVER_2D(m:m+l,'electrostatic_potential', savename) ;
       end
       disp([num2str(m) 'to' num2str(m+l) ' done'])
    end
end

%% 3. Process results

savename= 'examples/example1_ranges/example1_ranges_result.mat';

if ~exist(savename, 'file')
    for w = 1:length(M)
        m = M(w);
        for l = L
            load(['examples/precalculated/silver_' num2str(m) 'to' num2str(m+l) '.mat'],'eff_SILVER','ratio')
            min_eff_SILVER(w,l) = (min(eff_SILVER));
            min_eff_GR(w,l) = min(efficiency_range(gr2D,m:m+l,'electrostatic_potential'));
            prcnt_increase(w,l) = ((min_eff_SILVER(w,l)/min_eff_GR(w,l))-1)*100;
            SILVER_ratio(w,l) = ratio;

        end
    end
    save(savename, 'prcnt_increase','SILVER_ratio','min_eff_GR', 'min_eff_SILVER','M', 'L')
else
    warning('using pre-calculated results')
    load(savename);
end

%% 4. Vizualise result
tiny_ratios = [];
tiny_ratios_higher_order = [];
for N = 1
    tiny_ratios = cat(1,tiny_ratios,gen_gr2D(N,1));
    tiny_ratios = cat(1,1-tiny_ratios,gen_gr2D(N,1));
end

figure(1)
c_map = get(gca,'colororder');
plot(L,min_eff_SILVER(1,:), 'o-','Linewidth', 3, 'markersize', 5, 'Color',c_map(1,:))
hold on
plot(L,min_eff_GR(1,:), '--','Linewidth', 2, 'markersize', 5, 'Color',[c_map(1,:), 0.5])
plot(L,min_eff_SILVER(2,:),'o-','Linewidth', 3,'markersize', 5, 'Color',c_map(2,:))
plot(L,min_eff_GR(2,:), '--','Linewidth', 2, 'markersize', 5, 'Color',[c_map(2,:), 0.5])
plot(L,min_eff_SILVER(3,:),'o-', 'Linewidth', 3,'markersize', 5, 'Color',c_map(3,:))
plot(L,min_eff_GR(3,:), '--','Linewidth', 2, 'markersize', 5, 'Color',[c_map(3,:), 0.5])

set(gca,'FontSize',18)
set(gcf,'Position',[124 359 876 439])
axis([0,max(L)+1,ylim])
legend('SILVER, M = 4', 'GR, M = 4', 'SILVER, M = 16', 'GR, M = 16', 'SILVER, M = 32', 'GR, M = 32','location', 'northeast')
xlabel('Range width, R')
ylabel('Minimum efficiency, \eta')


savefig('examples/example1_ranges/example1_ranges_result.fig')
saveas(gcf,'examples/example1_ranges/example1_ranges_result.tiff')


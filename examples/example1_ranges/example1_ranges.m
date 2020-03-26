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
            SILVER_2D(m:m+l,'electrostatic_potential',savename) ;
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
            load(['examples/precalculated/silver_' num2str(m) 'to' num2str(m+l) '.mat'],'eff_SILVER')
            prcnt_increase(w,l) = (min(eff_SILVER)/min(efficiency_range(gr2D,m:m+l,'electrostatic_potential'))-1)*100;

        end
    end
    save(savename, 'prcnt_increase','M', 'L')
else
    warning('using pre-calculated results')
    load(savename);
end

%% 4. Vizualise result

figure(1)
plot(L+1,prcnt_increase(1,:), 'x-','Linewidth', 3, 'markersize', 10)
hold on
plot(L+1,prcnt_increase(2,:),'x-','Linewidth', 3,'markersize', 10)
plot(L+1,prcnt_increase(3,:),'x-', 'Linewidth', 3,'markersize', 10)

set(gca,'FontSize',18)
set(gcf,'Position',[124 359 876 439])
axis([0,max(L)+1,0,5])
legend('Minimum window size = 4', 'Minimum window size = 16', 'Minimum window size = 32', 'location', 'northeast')
xlabel('Range width')
ylabel('SILVER efficiency increase, %')

savefig('examples/example1_ranges/example1_ranges_result.fig')
saveas(gcf,'examples/example1_ranges/example1_ranges_result.tiff')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   How does efficiency, as calculated by the electrostatic potential, vary
%   within a range of window sizes or for multiple temporal resolutions for
%   SILVER and for the golden ratio method?
%                                              
%   Sophie Schauman 2020                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% 1. Choose  sets of window sizes, S, to consider

S1 = [16,32,48];
S2 = [33:39];
S = {S1,S2};

%% 2. Do SILVER optimization for those ranges
for s = S
   savename = ['examples/precalculated/silver_' strrep(num2str(s{:}),' ', '_') '.mat'];
   if ~exist(savename, 'file')
        SILVER_2D(s{:},'electrostatic_potential',savename) ;
   end
   disp([num2str(s{:}) ' done'])
end


%% 3. Vizualise results


eff_GR = efficiency_range(gr2D,10:50,'electrostatic_potential');

figure(3)
c_map = get(gca,'colororder');

hold on



n = 0;
for s = S
    n = n+1;
    for m = 1:length(s{:})
        rectangle('Position',[s{:}(m)-0.5,0.9,1,0.1], 'EdgeColor', 'None', 'FaceColor', [c_map(n,:), 0.5])
    end
    
    
    load(['examples/precalculated/silver_' strrep(num2str(s{:}),' ', '_') '.mat'], 'ratio')
    eff_SILVER_all = efficiency_range(ratio,10:50,'electrostatic_potential');
    h = plot(10:50,eff_SILVER_all,'x-','Linewidth', 3, 'markersize', 5);
    lg{n} = ['SILVER: S = \{' num2str(s{:}) '\}, \alpha = ' num2str(ratio)];

end
plot(10:50,eff_GR,'x-','Linewidth', 3, 'markersize', 5)
lg{end+1} = ['Golden ratio: \alpha = ' num2str(gr2D)];

box on
set(gca,'FontSize',18)
set(gcf,'Position',[124 359 876 439])
legend(lg{:}, 'location', 'southeast')
xlabel('Window size')
ylabel('Efficiency, \eta')
axis([xlim,0.9,1])

savefig('examples/example3_ep_efficiency/example3_ep_efficiency_result.fig')
saveas(gcf,'examples/example3_ep_efficiency/example3_ep_efficiency_result.tiff')


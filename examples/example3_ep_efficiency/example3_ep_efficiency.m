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

S1 = [15,25,50];
S2 = [32:39];
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


eff_GR = efficiency_range(gr2D,2:50,'electrostatic_potential');

figure(3)
plot(2:50,eff_GR,'-o','Linewidth', 3, 'markersize', 10)
hold on

lg{1} = 'Golden ratio';
n = 1;
for s = S
    n = n+1;
    load(['examples/precalculated/silver_' strrep(num2str(s{:}),' ', '_') '.mat'],'eff_SILVER', 'ratio')
    eff_SILVER_all = efficiency_range(ratio,2:50,'electrostatic_potential');
    h = plot(2:50,eff_SILVER_all,'x-','Linewidth', 3, 'markersize', 5);
    c = get(h,'Color');
    h.Color(4) = 0.2;
    plot(s{:},eff_SILVER,'o','Linewidth', 3, 'markersize', 10,'color', c)
    lg{n} = ['SILVER: S = \{' num2str(s{:}) '\} full range'];
    n = n+1;
    lg{n} = ['SILVER: S = \{' num2str(s{:}) '\} in optimized region'];
end

set(gca,'FontSize',18)
set(gcf,'Position',[124 359 876 439])
legend(lg{:}, 'location', 'southoutside')
xlabel('Window size')
ylabel('Efficiency')
axis([0,52,0.9,1])

savefig('examples/example3_ep_efficiency/example3_ep_efficiency_result.fig')
saveas(gcf,'examples/example3_ep_efficiency/example3_ep_efficiency_result.tiff')


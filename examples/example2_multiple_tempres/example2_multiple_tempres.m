%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Does SILVER perform better than the golden ratio when the set of window
%   sizes to optimize for are discrete temporal resolutions?
%                                              
%   Sophie Schauman 2020                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all


%% 1. Choose  sets of window sizes, S, to consider

S1 = [4,8];
S2 = [16,32];
S3 = [32,64];
S4 = [4,8,12];
S5 = [16,32,48];
S6 = [32,64,96];
S7 = fibonacci(5:9);
S = {S1,S2,S3,S4,S5,S6,S7};

%% 2. Do SILVER optimization for those ranges
for s = S
   savename = ['examples/precalculated/silver_' strrep(num2str(s{:}),' ', '_') '.mat'];
   if ~exist(savename, 'file')
        SILVER_2D(s{:},'electrostatic_potential',savename) ;
   end
   disp([num2str(s{:}) ' done'])
end


%% 3. Process results

savename= 'examples/example2_multiple_tempres/example2_multiple_tempres_result.mat';

if ~exist(savename, 'file')
    n = 0;
    for s = S
        n = n+1;
        load(['examples/precalculated/silver_' strrep(num2str(s{:}),' ', '_') '.mat'],'eff_SILVER','ratio')
            min_eff_SILVER(n) = (min(eff_SILVER));
            min_eff_GR(n) = min(efficiency_range(gr2D,s{:},'electrostatic_potential'));
            prcnt_increase(n) = ((min_eff_SILVER(n)/min_eff_GR(n))-1)*100;
            SILVER_ratio(n) = ratio;
    end
    save(savename, 'prcnt_increase','SILVER_ratio','min_eff_GR', 'min_eff_SILVER','S')
else
    warning('using pre-calculated results')
    load(savename);
end

%% 4. Vizualise result

figure(2)
c_map = [0.5 0.5 0.5; 1 0.5 0];

for n = 1:length(S)
    labels{n} = ['S = \{' num2str(S{n}) '\}'];
end
b = bar([min_eff_SILVER;min_eff_GR]','barwidth', 1);
b(1).FaceColor = c_map(1,:);
b(2).FaceColor = c_map(2,:);
b(1).LineWidth = 2;
b(2).LineWidth = 2;

set(gca,'xticklabel',labels)
xtickangle(-90)
legend('SILVER', 'Golden ratio','location', 'northeast')

axis([xlim, 0.9,1.01])
set(gca,'FontSize',16)
set(gca, 'LineWidth', 2)
grid on
set(gcf,'Position',[124 357 876 441])

ylabel('Minimum efficiency, \eta')

savefig('examples/example2_multiple_tempres/example2_multiple_tempres_result.fig')
saveas(gcf,'examples/example2_multiple_tempres/example2_multiple_tempres_result.tiff')

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
        load(['examples/precalculated/silver_' strrep(num2str(s{:}),' ', '_') '.mat'],'eff_SILVER')
        prcnt_increase(n) = (min(eff_SILVER)/min(efficiency_range(gr2D,s{:},'electrostatic_potential'))-1)*100;
    end
    save(savename, 'prcnt_increase','S')
else
    warning('using pre-calculated results')
    load(savename);
end

%% 4. Vizualise result

figure(2)
for n = 1:length(S)
    labels{n} = ['S = \{' num2str(S{n}) '\}'];
end
bar(prcnt_increase,'barwidth', 0.5)
set(gca,'xticklabel',labels)
xtickangle(-90)

set(gca,'FontSize',18)
set(gcf,'Position',[124 359 876 439])

ylabel('SILVER efficiency increase, %')

savefig('examples/example2_multiple_tempres/example2_multiple_tempres_result.fig')
saveas(gcf,'examples/example2_multiple_tempres/example2_multiple_tempres_result.tiff')
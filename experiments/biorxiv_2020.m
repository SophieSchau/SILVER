%% Figure 1
% Generated outside Matlab (not a results figure).

%% Figure 2

%%%%%%%% Parameters %%%%%%%%

S1 = [16,32,48];
S2 = [33:39];
S = {S1,S2};

xrange = 10:50;
ylims = [0.9,1];

%%%%%%%% Experiment + Figure generation %%%%%%%%
figure
SILVER_vs_GR_optimisation_example_fig(S, 'electrostatic_potential', xrange, ylims,'experiments/data/experiment_results/biorxiv2020/Figure2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Figure 3

%%%%%%%% Parameters %%%%%%%%

M = [4,16,32];
R = 1:100;

S1 = [4,8];
S2 = [16,32];
S3 = [32,64];
S4 = [4,8,12];
S5 = [16,32,48];
S6 = [32,64,96];
S7 = fibonacci(5:9);
S = {S1,S2,S3,S4,S5,S6,S7};

%%%%%%%% Experiment + Figure generation %%%%%%%%
close all
figure
SILVER_vs_GR_ranges(M, R, 'electrostatic_potential', 'experiments/data/experiment_results/biorxiv2020/Figure3a');
figure
SILVER_vs_GR_discrete_sets(S, 'electrostatic_potential', 'experiments/data/experiment_results/biorxiv2020/SILVERvsGR');
SILVER_vs_GR_discrete_sets_figure(S, 'experiments/data/experiment_results/biorxiv2020/SILVERvsGR', 'experiments/data/experiment_results/biorxiv2020/Figure3b');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Figure 4 

% NOTE: this method has some weaknesses that were kindly pointed out during 
% peer review that will be addressed prior to any further publications.

g_factors_biorxiv2020_script 



%% Figure 5

% NOTE: these methods have some weaknesses that were kindly pointed out during 
% peer review that will be addressed prior to any further publications.

phantom_biorxiv2020_script

clear
Subj = 1;
biorxiv2020_invivo_single_subj

clear
Subj = 2;
biorxiv2020_invivo_single_subj

clear
Subj = 3;
biorxiv2020_invivo_single_subj

biorxiv2020_invivo_group_analysis



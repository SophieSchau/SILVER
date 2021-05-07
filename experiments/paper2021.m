%% Figure 1
% Generated outside Matlab (not a results figure).

%% Figure 2
% Generated outside Matlab (not a results figure).

%% Figure 3a

%%%%%%%% Parameters %%%%%%%%

% discrete sets
S1 = [4,8];
S2 = [16,32];
S3 = [32,64];
S4 = [4,8,12];
S5 = [16,32,48];
S6 = [32,64,96];
S7 = fibonacci(5:9);

% continuous ranges
M = [8,16,32];
R = [3,5,7];

S8 = [M(1):M(1)+R(1)];
S9 = [M(1):M(1)+R(2)];
S10 = [M(1):M(1)+R(3)];

S11 = [M(2):M(2)+R(1)];
S12 = [M(2):M(2)+R(2)];
S13 = [M(2):M(2)+R(3)];

S14 = [M(3):M(3)+R(1)];
S15 = [M(3):M(3)+R(2)];
S16 = [M(3):M(3)+R(3)];

S_A = {S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16};

%%%%%%%% Experiment %%%%%%%%

SILVER_vs_GR_discrete_sets(S_A, 'electrostatic_potential', 'experiments/data/experiment_results/paper2021/SILVER_vs_GR_sets');

%%%%%%%% Figure %%%%%%%%
close all
order_sets = SILVER_vs_GR_discrete_sets_figure(S_A, 'experiments/data/experiment_results/paper2021/SILVER_vs_GR_sets', 'experiments/data/experiment_results/paper2021/figures/Figure3a');
S_A = S_A(order_sets);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 3b

%%%%%%%% Parameters %%%%%%%%

S1 = [16,32];
S2 = [16:19];
S_B = {S1,S2};

xrange = 10:40;
ylims = [0.9,1];

%%%%%%%% Figure %%%%%%%%
close all
SILVER_vs_GR_optimisation_example_fig(S_B, 'electrostatic_potential', xrange, ylims,'experiments/data/experiment_results/paper2021/figures/Figure3b');



%% Figure 3 - merge a and b
close all
figure
h(1)= subplot(2,1,1);
SILVER_vs_GR_discrete_sets_figure(S_A, 'experiments/data/experiment_results/paper2021/SILVER_vs_GR_sets')

%%%% Bold any sets highlighed in B

ax = gca;
for iA = 1:length(S_A)
    for iB = 1:length(S_B)
        if ismember(num2str(S_A{iA}), {num2str(S_B{iB})})
            ax.XTickLabel{iA} = ['\bf{' ax.XTickLabel{iA}, '}'];
        end
    end  
end

%%%%

h(2)=subplot(2,1,2);
SILVER_vs_GR_optimisation_example_fig(S_B, 'electrostatic_potential', xrange, ylims);

annotation('textbox',get(h(1), 'Position')*0.1+[0 0.9 0 0], 'String', '(A)', 'EdgeColor', 'none', 'fontsize', 16)
annotation('textbox',get(h(1), 'Position')*0.1+[0 0.4 0 0], 'String', '(B)', 'EdgeColor', 'none', 'fontsize', 16)

set(gcf, 'Position', [124 6 952 792])

savefig(gcf, 'experiments/data/experiment_results/paper2021/figures/Figure3ab')
saveas(gcf, 'experiments/data/experiment_results/paper2021/figures/Figure3ab.svg')


%% Figure 4
clear
close all
%%%%%%%% Parameters - single coil experiment %%%%%%%%
sens = ones(30);
ncov  = 1;
spokes = 60;
SILVER_sets = {[spokes-10, spokes], [spokes, spokes+10], [spokes-2, spokes, spokes+2]};
start_angle = 0; %degrees
%%%%%%%%% experiment - single coil

measure_noise_amplification(sens, ncov, spokes, SILVER_sets, 'electrostatic_potential', start_angle, ['experiments/data/experiment_results/paper2021/noise_single_coil_'])

%%%%%%%% figure
figure
measure_noise_amplification_single_coil_fig(spokes, ['experiments/data/experiment_results/paper2021/noise_single_coil*N_' num2str(spokes) '*.mat'], 'experiments/data/experiment_results/paper2021/figures/Figure4')



%% Figure 5
clear
close all
%%%%%%%% Parameters - with coils experiment %%%%%%%%
load('experiments/data/experiment_inputs/sensitivity_maps/meas_MID00130_FID121040_mc_ep3d_turbine_1a_sens0.3.mat','psens')
sens = psens(:,:,4,:); % use compressed coils
sensmask = logical(abs(sum(sens.^2,4)));
ncov  = 1;
spokes = 10;
SILVER_sets = {[5,10], [10,15], [9:11]};

%%%%%%%%% experiment - with coils

measure_noise_amplification(sens, ncov, spokes, SILVER_sets, 'electrostatic_potential', 0, ['experiments/data/experiment_results/paper2021/noise_multi_coil_sa_0_'])


%%%%%%%% figure
figure
measure_noise_amplification_multi_coil_fig(spokes, sensmask, ['experiments/data/experiment_results/paper2021/noise_multi_coil_sa_0*N_' num2str(spokes) '*.mat'], 'experiments/data/experiment_results/paper2021/figures/Figure5')

%% Figure 6
clear
close all
load('experiments/data/experiment_inputs/sensitivity_maps/meas_MID00130_FID121040_mc_ep3d_turbine_1a_sens0.3.mat','psens')
sens = psens(:,:,4,:); % use compressed coils
sensmask = logical(abs(sum(sens.^2,4)));
ncov  = 1;
spokes = 10;
SILVER_sets = {[5,10], [10,15], [9:11]};
start_angles = 1:179; %degrees
%%%%%%%%% experiment - with coils - multiple start angles

for start_angle = start_angles
    measure_noise_amplification(sens, ncov, spokes, SILVER_sets, 'electrostatic_potential', start_angle, ['experiments/data/experiment_results/paper2021/noise_multi_coil_sa_' num2str(start_angle) '_'])
end

%%%%%%%%%%% figure 

multi_start_angle_fig(start_angles, sensmask, 'experiments/data/experiment_results/paper2021/noise_multi_coil_sa_', 'experiments/data/experiment_results/paper2021/figures/Figure6')

%% Figure 7
% TURBINE acquisition in phantom - subtracted acquisitions
clear 
close all
%%%%%%%% SETUP files %%%%%%%% 
% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20201215_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/')
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20201215_TURBINE/meas_MID00130_FID121040_mc_ep3d_turbine_1a.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/meas_MID00130_FID121040_mc_ep3d_turbine_1a_slice',0.3, 2048,1:8);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/meas_MID00130_FID121040_mc_ep3d_turbine_1a_slice',1, 2048,1:8);
% recon_phantom(1:8,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/')


% insert option to data download here for when twix files are not available
load('experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/covariance.mat', 'ncov');
load('experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/meas_MID00130_FID121040_mc_ep3d_turbine_1a_sens0.3.mat','sens');
sens030 = sens;
load('experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/meas_MID00130_FID121040_mc_ep3d_turbine_1a_sens1.mat','sens');
sens100 = sens;

%%%%%%% Parameters
slices_quali = 4;
slices_quanti = 1:8;
sensmask_lowres = logical(abs(sum(sens030(:,:,:,:).^2,4)));
sensmask_highres = logical(abs(sum(sens100(:,:,:,:).^2,4)));
Nspokes = 2048;

%%%%%%% Expertiment
SNR_predict_and_measure_phantom(slices_quanti, Nspokes, sens030, sens100, ncov, 'experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/recons/')
SNR_predict_and_measure_phantom(slices_quali, Nspokes, sens030, sens100, ncov, 'experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/recons/')

%%%%%%% Figure
figure;
SNR_quali_fig(reshape(sensmask_lowres(:,:,slices_quali),30,[]), reshape(sensmask_highres(:,:,slices_quali),100,[]), slices_quali, 'experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/recons/', 'experiments/data/experiment_results/paper2021/figures/Figure7a')

figure;
SNR_quanti_fig(reshape(sensmask_lowres(:,:,slices_quanti),30,[]), reshape(sensmask_highres(:,:,slices_quanti),100,[]), slices_quanti, 'experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/recons/', 'experiments/data/experiment_results/paper2021/figures/Figure7b')

figure
A = imread('experiments/data/experiment_results/paper2021/figures/Figure7a.tiff');
B = imread('experiments/data/experiment_results/paper2021/figures/Figure7b.tiff');
imagesc(cat(2, A,B));
axis image 
axis off
set(gca, 'Position', [0 0 1 1])
set(gcf, 'Position', [1 7 1214 798])
annotation('textbox',[0.05 0.88 0.1 0.1], 'String', '(A)', 'LineStyle', 'none', 'fontsize', 16)
annotation('textbox',[0.3 0.88 0.1 0.1], 'String', '(B)', 'LineStyle', 'none','fontsize', 16)

saveas(gcf, 'experiments/data/experiment_results/paper2021/figures/Figure7ab.tiff')



%% Figure 8
% TURBINE acquisition in vivo - subtracted acquisitions
clear 
close all
%%%%%%%% SETUP files %%%%%%%% 
% insert option to pre-processed data download here for when twix files are not available

%%%%%%%%%%%%
% Subj A
% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20210310_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/')% 
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00690_FID123364_mc_ep3d_turbine_GRa_slice',0.3, 1536,1:16);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00690_FID123364_mc_ep3d_turbine_GRa_slice',1, 1536,1:16);
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20210310_TURBINE/meas_MID00690_FID123364_mc_ep3d_turbine_GRa.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));
% recon_invivo(1:16,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/')

% load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/covariance.mat', 'ncov');
% ncov = ncov;
% load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00690_FID123364_mc_ep3d_turbine_GRa__sens0.3.mat','sens');
% sens030 = sens;
% load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00690_FID123364_mc_ep3d_turbine_GRa__sens1.mat','sens');
% sens100 = sens;
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
% Subj B
% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20210316_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00995_FID123669_mc_ep3d_turbine_GRa_slice',0.3, 1536,1:16);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00995_FID123669_mc_ep3d_turbine_GRa_slice',1, 1536,1:16);
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20210316_TURBINE/meas_MID00995_FID123669_mc_ep3d_turbine_GRa.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));
% recon_invivo(1:16,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')

load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/covariance.mat', 'ncov');
ncov = ncov;
load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00995_FID123669_mc_ep3d_turbine_GRa__sens0.3.mat','sens');
sens030 = sens;
load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00995_FID123669_mc_ep3d_turbine_GRa__sens1.mat','sens');
sens100 = sens;
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% Parameters
slices_quali = 4;
slices_quanti = 1:16;
Nspokes = 1536;
sensmask_lowres = logical(abs(sum(sens030(:,:,:,:).^2,4)));
sensmask_highres = logical(abs(sum(sens100(:,:,:,:).^2,4)));
subj = 'B';

%%%%%%% Expertiment
SNR_predict_and_measure_invivo(slices_quanti, Nspokes, sens030, sens100, ncov, ['experiments/data/experiment_inputs/TURBINE_data/in-vivo/subj' subj '/recons/'])
SNR_predict_and_measure_invivo(slices_quali, Nspokes, sens030, sens100, ncov, ['experiments/data/experiment_inputs/TURBINE_data/in-vivo/subj' subj '/recons/'])

%%%%%%% Figure
figure %A
SNR_quali_fig(reshape(sensmask_lowres(:,:,slices_quali),30,[]), reshape(sensmask_highres(:,:,slices_quali),100,[]), slices_quali, ['experiments/data/experiment_inputs/TURBINE_data/in-vivo/subj' subj '/recons/'], 'experiments/data/experiment_results/paper2021/figures/Figure8a')

figure %B
SNR_quanti_fig(reshape(sensmask_lowres(:,:,slices_quanti),30,[]), reshape(sensmask_highres(:,:,slices_quanti),100,[]), slices_quanti, ['experiments/data/experiment_inputs/TURBINE_data/in-vivo/subj' subj '/recons/'], 'experiments/data/experiment_results/paper2021/figures/Figure8b')

figure % Combine A and B
A = imread('experiments/data/experiment_results/paper2021/figures/Figure8a.tiff');
B = imread('experiments/data/experiment_results/paper2021/figures/Figure8b.tiff');
imagesc(cat(2, A,B));
axis image 
axis off
set(gca, 'Position', [0 0 1 1])
set(gcf, 'Position', [1 7 1214 798])
annotation('textbox',[0.05 0.88 0.1 0.1], 'String', '(A)', 'LineStyle', 'none', 'fontsize', 16)
annotation('textbox',[0.3 0.88 0.1 0.1], 'String', '(B)', 'LineStyle', 'none','fontsize', 16)

saveas(gcf, 'experiments/data/experiment_results/paper2021/figures/Figure8ab.tiff')


%% Figure 9
% TURBINE acquisition in vivo - resting state linear recon
clear 
close all
%%%%%%%% SETUP files %%%%%%%% 
% insert option to pre-processed data download here for when twix files are not available

%%%%%%%%%%%%
% Subj A
% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20210310_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/')% 
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI_slice',0.3, 6016,1:16);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI_slice',1, 6016,1:16);
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20210310_TURBINE/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));
% recon_invivo_fmri(1:16,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/')

load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/covariance_rs_fMRI.mat', 'ncov');
ncov = ncov;
load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI__sens0.3.mat','sens');
sens030 = sens;
load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI__sens1.mat','sens');
sens100 = sens;
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
% Subj B
% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20210316_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')% 
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI_slice',0.3, 6016,1:16);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI_slice',1, 6016,1:16);
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20210316_TURBINE/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));
% recon_invivo_fmri(1:16,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')

% load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/covariance_rs_fMRI.mat', 'ncov');
% ncov = ncov;
% load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI__sens0.3.mat','sens');
% sens030 = sens;
% load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI__sens1.mat','sens');
% sens100 = sens;
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% Parameters
slices_quali = 4;
slices_quanti = 3:13;
sensmask_lowres = logical(abs(sum(sens030(:,:,:,:).^2,4)));
sensmask_highres = logical(abs(sum(sens100(:,:,:,:).^2,4)));
subj = 'A';

%%%%%%% Figure
tSNR_fig(['experiments/data/experiment_inputs/TURBINE_data/in-vivo/subj' subj '/recons/'],3:13, 8, sensmask_lowres, sensmask_highres, 'experiments/data/experiment_results/paper2021/figures/Figure9')

%% Figure 10
% TURBINE acquisition in vivo - resting state nonlinear
clear 
close all
%%%%%%%% SETUP files %%%%%%%% 
% insert option to pre-processed data download here for when twix files are not available

%%%%%%%%%%%%
% Subj A
% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20210310_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/')% 
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI_slice',0.3, 6016,1:16);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI_slice',1, 6016,1:16);
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20210310_TURBINE/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));
% recon_invivo_fmri_wavelet(1:16,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/')

load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/covariance_rs_fMRI.mat', 'ncov');
ncov = ncov;
load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI__sens0.3.mat','sens');
sens030 = sens;
load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI__sens1.mat','sens');
sens100 = sens;
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
% Subj B
% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20210316_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')% 
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI_slice',0.3, 6016,1:16);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI_slice',1, 6016,1:16);
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20210316_TURBINE/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));
% recon_invivo_fmri_wavelet(1:16,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')

% load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/covariance_rs_fMRI.mat', 'ncov');
% ncov = ncov;
% load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI__sens0.3.mat','sens');
% sens030 = sens;
% load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI__sens1.mat','sens');
% sens100 = sens;
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% Parameters
slices_quali = 4;
slices_quanti = 3:13;
sensmask_lowres = logical(abs(sum(sens030(:,:,:,:).^2,4)));
sensmask_highres = logical(abs(sum(sens100(:,:,:,:).^2,4)));
subj = 'A';

%%%%%%% Figure
tSNR_fig(['experiments/data/experiment_inputs/TURBINE_data/in-vivo/subj' subj '/recons_nl_wavelet/'],3:13, 8, sensmask_lowres, sensmask_highres, 'experiments/data/experiment_results/paper2021/figures/Figure10')


%% Supplementary figure 1
%%%%%%%% Parameters %%%%%%%%
NSpokes = 60;
mat_size  = 30;
rng(1)
alphas = [gr2D, 1/NSpokes, rand(1, 100)*0.9+0.05];
%%%%%%%% Experiment + Figure generation %%%%%%%%
[efficiency,noise,alpha_sorted] = ...
efficiency_and_noise_amplification(NSpokes,...
                                   'electrostatic_potential',...
                                   mat_size,...
                                   ['experiments/data/experiment_results/paper2021/eff_and_noise_' num2str(NSpokes) '_' num2str(mat_size)],...
                                   alphas);
                               
savefig(['experiments/data/experiment_results/paper2021/figures/FigureS1.fig'])
saveas(gcf,['experiments/data/experiment_results/paper2021/figures/FigureS1.svg'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Supplementary figure 2
% TURBINE acquisition in vivo - subtracted acquisitions
clear 
close all
%%%%%%%% SETUP files %%%%%%%% 
% insert option to pre-processed data download here for when twix files are not available

%%%%%%%%%%%%
% Subj A
% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20210310_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/')% 
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00690_FID123364_mc_ep3d_turbine_GRa_slice',0.3, 1536,1:16);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00690_FID123364_mc_ep3d_turbine_GRa_slice',1, 1536,1:16);
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20210310_TURBINE/meas_MID00690_FID123364_mc_ep3d_turbine_GRa.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));
% recon_invivo(1:16,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/')

load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/covariance.mat', 'ncov');
ncov = ncov;
load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00690_FID123364_mc_ep3d_turbine_GRa__sens0.3.mat','sens');
sens030 = sens;
load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00690_FID123364_mc_ep3d_turbine_GRa__sens1.mat','sens');
sens100 = sens;
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
% Subj B
% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20210316_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00995_FID123669_mc_ep3d_turbine_GRa_slice',0.3, 1536,1:16);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00995_FID123669_mc_ep3d_turbine_GRa_slice',1, 1536,1:16);
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20210316_TURBINE/meas_MID00995_FID123669_mc_ep3d_turbine_GRa.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));
% recon_invivo(1:16,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')

% load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/covariance.mat', 'ncov');
% ncov = ncov;
% load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00995_FID123669_mc_ep3d_turbine_GRa__sens0.3.mat','sens');
% sens030 = sens;
% load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00995_FID123669_mc_ep3d_turbine_GRa__sens1.mat','sens');
% sens100 = sens;
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% Parameters
slices_quali = 4;
slices_quanti = 1:16;
Nspokes = 1536;
sensmask_lowres = logical(abs(sum(sens030(:,:,:,:).^2,4)));
sensmask_highres = logical(abs(sum(sens100(:,:,:,:).^2,4)));
subj = 'A';

%%%%%%% Expertiment
SNR_predict_and_measure_invivo(slices_quanti, Nspokes, sens030, sens100, ncov, ['experiments/data/experiment_inputs/TURBINE_data/in-vivo/subj' subj '/recons/'])
SNR_predict_and_measure_invivo(slices_quali, Nspokes, sens030, sens100, ncov, ['experiments/data/experiment_inputs/TURBINE_data/in-vivo/subj' subj '/recons/'])

%%%%%%% Figure
%%%%%%% Figure
figure %A
SNR_quali_fig(reshape(sensmask_lowres(:,:,slices_quali),30,[]), reshape(sensmask_highres(:,:,slices_quali),100,[]), slices_quali, ['experiments/data/experiment_inputs/TURBINE_data/in-vivo/subj' subj '/recons/'], 'experiments/data/experiment_results/paper2021/figures/Figure8a')

figure %B
SNR_quanti_fig(reshape(sensmask_lowres(:,:,slices_quanti),30,[]), reshape(sensmask_highres(:,:,slices_quanti),100,[]), slices_quanti, ['experiments/data/experiment_inputs/TURBINE_data/in-vivo/subj' subj '/recons/'], 'experiments/data/experiment_results/paper2021/figures/Figure8b')

figure % Combine A and B
A = imread('experiments/data/experiment_results/paper2021/figures/FigureS2a.tiff');
B = imread('experiments/data/experiment_results/paper2021/figures/FigureS2b.tiff');
imagesc(cat(2, A,B));
axis image 
axis off
set(gca, 'Position', [0 0 1 1])
set(gcf, 'Position', [1 7 1214 798])
annotation('textbox',[0.05 0.88 0.1 0.1], 'String', '(A)', 'LineStyle', 'none', 'fontsize', 16)
annotation('textbox',[0.3 0.88 0.1 0.1], 'String', '(B)', 'LineStyle', 'none','fontsize', 16)

saveas(gcf, 'experiments/data/experiment_results/paper2021/figures/FigureS2ab.tiff')
                               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary figure 3
% TURBINE acquisition in vivo - tSNR resting state
clear 
close all
%%%%%%%% SETUP files %%%%%%%% 
% insert option to pre-processed data download here for when twix files are not available

%%%%%%%%%%%%
% Subj A
% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20210310_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/')% 
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI_slice',0.3, 6016,1:16);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI_slice',1, 6016,1:16);
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20210310_TURBINE/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));
% recon_invivo_fmri(1:16,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/')

% load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/covariance_rs_fMRI.mat', 'ncov');
% ncov = ncov;
% load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI__sens0.3.mat','sens');
% sens030 = sens;
% load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI__sens1.mat','sens');
% sens100 = sens;
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
% Subj B
% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20210316_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')% 
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI_slice',0.3, 6016,1:16);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI_slice',1, 6016,1:16);
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20210316_TURBINE/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));
% recon_invivo_fmri(1:16,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')

load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/covariance_rs_fMRI.mat', 'ncov');
ncov = ncov;
load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI__sens0.3.mat','sens');
sens030 = sens;
load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI__sens1.mat','sens');
sens100 = sens;
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% Parameters
slices_quali = 4;
slices_quanti = 3:13;
sensmask_lowres = logical(abs(sum(sens030(:,:,:,:).^2,4)));
sensmask_highres = logical(abs(sum(sens100(:,:,:,:).^2,4)));
subj = 'B';

%%%%%%% Figure
tSNR_fig(['experiments/data/experiment_inputs/TURBINE_data/in-vivo/subj' subj '/recons/'],3:13, 8, sensmask_lowres, sensmask_highres, 'experiments/data/experiment_results/paper2021/figures/FigureS3')

%% Supplementary figure 4
% TURBINE acquisition in vivo - tSNR resting state nonlinear
clear 
close all
%%%%%%%% SETUP files %%%%%%%% 
% insert option to pre-processed data download here for when twix files are not available

%%%%%%%%%%%%
% Subj A
% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20210310_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/')% 
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI_slice',0.3, 6016,1:16);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI_slice',1, 6016,1:16);
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20210310_TURBINE/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));
% recon_invivo_fmri_wavelet(1:16,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/')

% load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/covariance_rs_fMRI.mat', 'ncov');
% ncov = ncov;
% load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI__sens0.3.mat','sens');
% sens030 = sens;
% load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI__sens1.mat','sens');
% sens100 = sens;
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
% Subj B
% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20210316_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')% 
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI_slice',0.3, 6016,1:16);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI_slice',1, 6016,1:16);
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20210316_TURBINE/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));
% recon_invivo_fmri_wavelet(1:16,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')

load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/covariance_rs_fMRI.mat', 'ncov');
ncov = ncov;
load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI__sens0.3.mat','sens');
sens030 = sens;
load('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI__sens1.mat','sens');
sens100 = sens;
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% Parameters
slices_quali = 4;
slices_quanti = 3:13;
sensmask_lowres = logical(abs(sum(sens030(:,:,:,:).^2,4)));
sensmask_highres = logical(abs(sum(sens100(:,:,:,:).^2,4)));
subj = 'B';

%%%%%%% Figure
tSNR_fig(['experiments/data/experiment_inputs/TURBINE_data/in-vivo/subj' subj '/recons_nl_wavelet/'],3:13, 8, sensmask_lowres, sensmask_highres, 'experiments/data/experiment_results/paper2021/figures/FigureS4')

%% Supplementary figure 5-8
% in-vivo reconstructions
%%%%%%% Parameters
clear
close all
%%%%% Parameters 5
slice = 8;
subj = 'A';
recon_folder =['experiments/data/experiment_inputs/TURBINE_data/in-vivo/subj' subj '/recons/'];
titl = ['Linear - Subj ' subj];
%%%%% Figure 5
figure
rs_recon_fig(recon_folder, slice, titl, 'experiments/data/experiment_results/paper2021/figures/FigureS5')


%%%%% Parameters 6
slice = 8;
subj = 'B';
titl = ['Linear - Subj ' subj];
recon_folder =['experiments/data/experiment_inputs/TURBINE_data/in-vivo/subj' subj '/recons/'];

%%%%% Figure
figure
rs_recon_fig(recon_folder, slice, titl, 'experiments/data/experiment_results/paper2021/figures/FigureS6')

%%%%% Parameters 7
slice = 8;
subj = 'A';
titl = ['Non-linear - Subj ' subj];
recon_folder =['experiments/data/experiment_inputs/TURBINE_data/in-vivo/subj' subj '/recons_nl_wavelet/'];

%%%%% Figure
figure
rs_recon_fig(recon_folder, slice, titl, 'experiments/data/experiment_results/paper2021/figures/FigureS7')

%%%%% Parameters 6
slice = 8;
subj = 'B';
titl = ['Non-linear - Subj ' subj];
recon_folder =['experiments/data/experiment_inputs/TURBINE_data/in-vivo/subj' subj '/recons_nl_wavelet/'];

%%%%% Figure
figure
rs_recon_fig(recon_folder, slice, titl, 'experiments/data/experiment_results/paper2021/figures/FigureS8')

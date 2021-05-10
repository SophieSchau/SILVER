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

if exist('experiments/data/experiment_inputs/sensitivity_maps/meas_MID00130_FID121040_mc_ep3d_turbine_1a_sens0.3.mat', 'file')
    load('experiments/data/experiment_inputs/sensitivity_maps/meas_MID00130_FID121040_mc_ep3d_turbine_1a_sens0.3.mat','psens')
else
    load(websave('experiments/data/experiment_inputs/sensitivity_maps/meas_MID00130_FID121040_mc_ep3d_turbine_1a_sens0.3.mat', 'https://zenodo.org/record/4743420/files/meas_MID00130_FID121040_mc_ep3d_turbine_1a_sens0.3.mat'),'psens')
end

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

if exist('experiments/data/experiment_inputs/sensitivity_maps/meas_MID00130_FID121040_mc_ep3d_turbine_1a_sens0.3.mat', 'file')
    load('experiments/data/experiment_inputs/sensitivity_maps/meas_MID00130_FID121040_mc_ep3d_turbine_1a_sens0.3.mat','psens')
else
    load(websave('experiments/data/experiment_inputs/sensitivity_maps/meas_MID00130_FID121040_mc_ep3d_turbine_1a_sens0.3.mat', 'https://zenodo.org/record/4743420/files/meas_MID00130_FID121040_mc_ep3d_turbine_1a_sens0.3.mat'),'psens')
end

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
% For reference of how the data was pre-processed:

% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20201215_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/')
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20201215_TURBINE/meas_MID00130_FID121040_mc_ep3d_turbine_1a.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/meas_MID00130_FID121040_mc_ep3d_turbine_1a_slice',0.3, 2048,1:8);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/meas_MID00130_FID121040_mc_ep3d_turbine_1a_slice',1, 2048,1:8);



% Data that is not available will be downloaded from Zenodo:

covfile_name = 'experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/covariance.mat';
sens_lowres_name = 'experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/meas_MID00130_FID121040_mc_ep3d_turbine_1a_sens0.3.mat';
sens_highres_name = 'experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/meas_MID00130_FID121040_mc_ep3d_turbine_1a_sens1.mat';

if exist(covfile_name, 'file')
    load(covfile_name, 'ncov');
else
    load(websave(covfile_name, 'https://zenodo.org/record/4743420/files/covariance.mat'), 'ncov');
end
    
if exist(sens_lowres_name, 'file')
    load(sens_lowres_name,'sens');
else
    load(websave(sens_lowres_name, 'https://zenodo.org/record/4743420/files/meas_MID00130_FID121040_mc_ep3d_turbine_1a_sens0.3.mat'), 'sens');
end
sens030 = sens;

if exist(sens_highres_name, 'file')
    load(sens_highres_name,'sens');
else
    load(websave(sens_highres_name, 'https://zenodo.org/record/4743420/files/meas_MID00130_FID121040_mc_ep3d_turbine_1a_sens1.mat'), 'sens');
end
sens100 = sens;


localpath = 'experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/';
webpath = 'https://zenodo.org/record/4743420/files/';

a1_name = 'meas_MID00130_FID121040_mc_ep3d_turbine_1a_slice';
b1_name = 'meas_MID00132_FID121042_mc_ep3d_turbine_1b_slice';
a2_name = 'meas_MID00134_FID121044_mc_ep3d_turbine_2a_slice';
b2_name = 'meas_MID00136_FID121046_mc_ep3d_turbine_2b_slice';
a3_name = 'meas_MID00138_FID121048_mc_ep3d_turbine_3a_slice';
b3_name = 'meas_MID00140_FID121050_mc_ep3d_turbine_3b_slice';
a4_name = 'meas_MID00142_FID121052_mc_ep3d_turbine_4a_slice';
b4_name = 'meas_MID00144_FID121054_mc_ep3d_turbine_4b_slice';
a5_name = 'meas_MID00146_FID121056_mc_ep3d_turbine_5a_slice';
b5_name = 'meas_MID00148_FID121058_mc_ep3d_turbine_5b_slice';
a6_name = 'meas_MID00150_FID121060_mc_ep3d_turbine_6a_slice';
b6_name = 'meas_MID00152_FID121062_mc_ep3d_turbine_6b_slice';
a8_name = 'meas_MID00158_FID121068_mc_ep3d_turbine_8a_slice';
b8_name = 'meas_MID00160_FID121070_mc_ep3d_turbine_8b_slice';

% download data
for slice = 1:8
    if ~exist([localpath a1_name num2str(slice) '.mat'], 'file')
        websave([localpath a1_name num2str(slice) '.mat'], [webpath a1_name num2str(slice) '.mat']);
    end
    if ~exist([localpath b1_name num2str(slice) '.mat'], 'file')
        websave([localpath b1_name num2str(slice) '.mat'], [webpath b1_name num2str(slice) '.mat']);
    end
    if ~exist([localpath a2_name num2str(slice) '.mat'], 'file')
        websave([localpath a2_name num2str(slice) '.mat'], [webpath a2_name num2str(slice) '.mat']);
    end
    if ~exist([localpath b2_name num2str(slice) '.mat'], 'file')
        websave([localpath b2_name num2str(slice) '.mat'], [webpath b2_name num2str(slice) '.mat']);
    end
    if ~exist([localpath a3_name num2str(slice) '.mat'], 'file')
        websave([localpath a3_name num2str(slice) '.mat'], [webpath a3_name num2str(slice) '.mat']);
    end
    if ~exist([localpath b3_name num2str(slice) '.mat'], 'file')
        websave([localpath b3_name num2str(slice) '.mat'], [webpath b3_name num2str(slice) '.mat']);
    end
    if ~exist([localpath a4_name num2str(slice) '.mat'], 'file')
        websave([localpath a4_name num2str(slice) '.mat'], [webpath a4_name num2str(slice) '.mat']);
    end
    if ~exist([localpath b4_name num2str(slice) '.mat'], 'file')
        websave([localpath b4_name num2str(slice) '.mat'], [webpath b4_name num2str(slice) '.mat']);
    end
    if ~exist([localpath a5_name num2str(slice) '.mat'], 'file')
        websave([localpath a5_name num2str(slice) '.mat'], [webpath a5_name num2str(slice) '.mat']);
    end
    if ~exist([localpath b5_name num2str(slice) '.mat'], 'file')
        websave([localpath b5_name num2str(slice) '.mat'], [webpath b5_name num2str(slice) '.mat']);
    end
    if ~exist([localpath a6_name num2str(slice) '.mat'], 'file')
        websave([localpath a6_name num2str(slice) '.mat'], [webpath a6_name num2str(slice) '.mat']);
    end
    if ~exist([localpath b6_name num2str(slice) '.mat'], 'file')
        websave([localpath b6_name num2str(slice) '.mat'], [webpath b6_name num2str(slice) '.mat']);
    end
    if ~exist([localpath a8_name num2str(slice) '.mat'], 'file')
        websave([localpath a8_name num2str(slice) '.mat'], [webpath a8_name num2str(slice) '.mat']);
    end
    if ~exist([localpath b8_name num2str(slice) '.mat'], 'file')
        websave([localpath b8_name num2str(slice) '.mat'], [webpath b8_name num2str(slice) '.mat']);
    end
end

% perform reconstruction (skips if files already exist)
recon_phantom(1:8,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/phantom/DEC15/')


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
% TURBINE acquisition in vivo - subtracted acquisitions subj B
clear 
close all
%%%%%%%% SETUP files %%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%
% Subj B - for reference preprocessing pipeline:
% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20210316_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00995_FID123669_mc_ep3d_turbine_GRa_slice',0.3, 1536,1:16);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00995_FID123669_mc_ep3d_turbine_GRa_slice',1, 1536,1:16);
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20210316_TURBINE/meas_MID00995_FID123669_mc_ep3d_turbine_GRa.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));
% recon_invivo(1:16,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')



% Data that is not available will be downloaded from Zenodo:

covfile_name = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/covariance.mat';
sens_lowres_name = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00995_FID123669_mc_ep3d_turbine_GRa__sens0.3.mat';
sens_highres_name = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00995_FID123669_mc_ep3d_turbine_GRa__sens1.mat';

if exist(covfile_name, 'file')
    load(covfile_name, 'ncov');
else
    load(websave(covfile_name, 'https://zenodo.org/record/4743764/files/covariance.mat'), 'ncov');
end
    
if exist(sens_lowres_name, 'file')
    load(sens_lowres_name,'sens');
else
    load(websave(sens_lowres_name, 'https://zenodo.org/record/4743764/files/meas_MID00995_FID123669_mc_ep3d_turbine_GRa__sens0.3.mat'), 'sens');
end
sens030 = sens;

if exist(sens_highres_name, 'file')
    load(sens_highres_name,'sens');
else
    load(websave(sens_highres_name, 'https://zenodo.org/record/4743764/files/meas_MID00995_FID123669_mc_ep3d_turbine_GRa__sens1.mat'), 'sens');
end
sens100 = sens;


localpath = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/';
webpath = 'https://zenodo.org/record/4743764/files/';

GRa_name = 'meas_MID00995_FID123669_mc_ep3d_turbine_GRa_slice';
GRb_name = 'meas_MID00997_FID123671_mc_ep3d_turbine_GRb_slice';
U8a_name = 'meas_MID01007_FID123681_mc_ep3d_turbine_UNI_8a_slice';
U8b_name = 'meas_MID01009_FID123683_mc_ep3d_turbine_UNI_8b_slice';
U10a_name = 'meas_MID00999_FID123673_mc_ep3d_turbine_UNI_10a_slice';
U10b_name = 'meas_MID01001_FID123675_mc_ep3d_turbine_UNI_10b_slice';
U46a_name = 'meas_MID01019_FID123693_mc_ep3d_turbine_UNI_46a_slice';
U46b_name = 'meas_MID01021_FID123695_mc_ep3d_turbine_UNI_46b_slice';
U55a_name = 'meas_MID01011_FID123685_mc_ep3d_turbine_UNI_55a_slice';
U55b_name = 'meas_MID01013_FID123687_mc_ep3d_turbine_UNI_55b_slice';
S10_46a_name = 'meas_MID01003_FID123677_mc_ep3d_turbine_SILVER_10_46a_slice';
S10_46b_name = 'meas_MID01005_FID123679_mc_ep3d_turbine_SILVER_10_46b_slice';
S8_55a_name = 'meas_MID01015_FID123689_mc_ep3d_turbine_SILVER_8_55a_slice';
S8_55b_name = 'meas_MID01017_FID123691_mc_ep3d_turbine_SILVER_8_55b_slice';

% download data
for slice = 1:16
    if ~exist([localpath GRa_name num2str(slice) '.mat'], 'file')
        websave([localpath GRa_name num2str(slice) '.mat'], [webpath GRa_name num2str(slice) '.mat']);
    end
    if ~exist([localpath GRb_name num2str(slice) '.mat'], 'file')
        websave([localpath GRb_name num2str(slice) '.mat'], [webpath GRb_name num2str(slice) '.mat']);
    end
    if ~exist([localpath U8a_name num2str(slice) '.mat'], 'file')
        websave([localpath U8a_name num2str(slice) '.mat'], [webpath U8a_name num2str(slice) '.mat']);
    end
    if ~exist([localpath U8b_name num2str(slice) '.mat'], 'file')
        websave([localpath U8b_name num2str(slice) '.mat'], [webpath U8b_name num2str(slice) '.mat']);
    end
    if ~exist([localpath U10a_name num2str(slice) '.mat'], 'file')
        websave([localpath U10a_name num2str(slice) '.mat'], [webpath U10a_name num2str(slice) '.mat']);
    end
    if ~exist([localpath U10b_name num2str(slice) '.mat'], 'file')
        websave([localpath U10b_name num2str(slice) '.mat'], [webpath U10b_name num2str(slice) '.mat']);
    end
    if ~exist([localpath U46a_name num2str(slice) '.mat'], 'file')
        websave([localpath U46a_name num2str(slice) '.mat'], [webpath U46a_name num2str(slice) '.mat']);
    end
    if ~exist([localpath U46b_name num2str(slice) '.mat'], 'file')
        websave([localpath U46b_name num2str(slice) '.mat'], [webpath U46b_name num2str(slice) '.mat']);
    end
    if ~exist([localpath U55a_name num2str(slice) '.mat'], 'file')
        websave([localpath U55a_name num2str(slice) '.mat'], [webpath U55a_name num2str(slice) '.mat']);
    end
    if ~exist([localpath U55b_name num2str(slice) '.mat'], 'file')
        websave([localpath U55b_name num2str(slice) '.mat'], [webpath U55b_name num2str(slice) '.mat']);
    end
    if ~exist([localpath S8_55a_name num2str(slice) '.mat'], 'file')
        websave([localpath S8_55a_name num2str(slice) '.mat'], [webpath S8_55a_name num2str(slice) '.mat']);
    end
    if ~exist([localpath S8_55b_name num2str(slice) '.mat'], 'file')
        websave([localpath S8_55b_name num2str(slice) '.mat'], [webpath S8_55b_name num2str(slice) '.mat']);
    end
    if ~exist([localpath S10_46a_name num2str(slice) '.mat'], 'file')
        websave([localpath S10_46a_name num2str(slice) '.mat'], [webpath S10_46a_name num2str(slice) '.mat']);
    end
    if ~exist([localpath S10_46b_name num2str(slice) '.mat'], 'file')
        websave([localpath S10_46b_name num2str(slice) '.mat'], [webpath S10_46b_name num2str(slice) '.mat']);
    end
end

% perform reconstruction (skips if files already exist)
recon_invivo(1:16,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')



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
% TURBINE acquisition in vivo - resting state linear recon subj A
clear 
close all
%%%%%%%% SETUP files %%%%%%%% 
%%%%%%%%%%%%
% Subj A
% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20210310_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/')% 
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI_slice',0.3, 6016,1:16);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI_slice',1, 6016,1:16);
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20210310_TURBINE/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));

covfile_name = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/covariance_rs_fMRI.mat';
sens_lowres_name = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI__sens0.3.mat';
sens_highres_name = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI__sens1.mat';

if exist(covfile_name, 'file')
    load(covfile_name, 'ncov');
else
    load(websave(covfile_name, 'https://zenodo.org/record/4743418/files/covariance.mat'), 'ncov');
end
    
if exist(sens_lowres_name, 'file')
    load(sens_lowres_name,'sens');
else
    load(websave(sens_lowres_name, 'https://zenodo.org/record/4743418/files/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI__sens0.3.mat'), 'sens');
end
sens030 = sens;

if exist(sens_highres_name, 'file')
    load(sens_highres_name,'sens');
else
    load(websave(sens_highres_name, 'https://zenodo.org/record/4743418/files/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI__sens1.mat'), 'sens');
end
sens100 = sens;


localpath = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/';
webpath = 'https://zenodo.org/record/4743418/files/';

GR_name = 'meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI_slice';
S_name = 'meas_MID00669_FID123343_mc_ep3d_turbine_SILVER_rs_fMRI_slice';

% download data
for slice = 3:13
    if ~exist([localpath GR_name num2str(slice) '.mat'], 'file')
        websave([localpath GR_name num2str(slice) '.mat'], [webpath GR_name num2str(slice) '.mat']);
    end
    if ~exist([localpath S_name num2str(slice) '.mat'], 'file')
        websave([localpath S_name num2str(slice) '.mat'], [webpath S_name num2str(slice) '.mat']);
    end
end

% reconstruct
recon_invivo_fmri(3:13,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/')


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

covfile_name = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/covariance_rs_fMRI.mat';
sens_lowres_name = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI__sens0.3.mat';
sens_highres_name = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI__sens1.mat';

if exist(covfile_name, 'file')
    load(covfile_name, 'ncov');
else
    load(websave(covfile_name, 'https://zenodo.org/record/4743418/files/covariance.mat'), 'ncov');
end
    
if exist(sens_lowres_name, 'file')
    load(sens_lowres_name,'sens');
else
    load(websave(sens_lowres_name, 'https://zenodo.org/record/4743418/files/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI__sens0.3.mat'), 'sens');
end
sens030 = sens;

if exist(sens_highres_name, 'file')
    load(sens_highres_name,'sens');
else
    load(websave(sens_highres_name, 'https://zenodo.org/record/4743418/files/meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI__sens1.mat'), 'sens');
end
sens100 = sens;


localpath = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/';
webpath = 'https://zenodo.org/record/4743418/files/';

GR_name = 'meas_MID00667_FID123341_mc_ep3d_turbine_GR_rs_fMRI_slice';
S_name = 'meas_MID00669_FID123343_mc_ep3d_turbine_SILVER_rs_fMRI_slice';

% download data
for slice = 3:13
    if ~exist([localpath GR_name num2str(slice) '.mat'], 'file')
        websave([localpath GR_name num2str(slice) '.mat'], [webpath GR_name num2str(slice) '.mat']);
    end
    if ~exist([localpath S_name num2str(slice) '.mat'], 'file')
        websave([localpath S_name num2str(slice) '.mat'], [webpath S_name num2str(slice) '.mat']);
    end
end

% reconstruct
recon_invivo_fmri_wavelet(1:16,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/')



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
% TURBINE acquisition in vivo - subtracted acquisitions subj A
clear 
close all
%%%%%%%% SETUP files %%%%%%%% 
% insert option to pre-processed data download here for when twix files are not available

%%%%%%%%%%%%
% Subj A
% For reference of how the data was pre-processed:

% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20210310_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/')% 
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00690_FID123364_mc_ep3d_turbine_GRa_slice',0.3, 1536,1:16);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00690_FID123364_mc_ep3d_turbine_GRa_slice',1, 1536,1:16);
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20210310_TURBINE/meas_MID00690_FID123364_mc_ep3d_turbine_GRa.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));



% Data that is not available will be downloaded from Zenodo:

covfile_name = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/covariance.mat';
sens_lowres_name = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00690_FID123364_mc_ep3d_turbine_GRa__sens0.3.mat';
sens_highres_name = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/meas_MID00690_FID123364_mc_ep3d_turbine_GRa__sens1.mat';

if exist(covfile_name, 'file')
    load(covfile_name, 'ncov');
else
    load(websave(covfile_name, 'https://zenodo.org/record/4743418/files/covariance.mat'), 'ncov');
end
    
if exist(sens_lowres_name, 'file')
    load(sens_lowres_name,'sens');
else
    load(websave(sens_lowres_name, 'https://zenodo.org/record/4743418/files/meas_MID00690_FID123364_mc_ep3d_turbine_GRa__sens0.3.mat'), 'sens');
end
sens030 = sens;

if exist(sens_highres_name, 'file')
    load(sens_highres_name,'sens');
else
    load(websave(sens_highres_name, 'https://zenodo.org/record/4743418/files/meas_MID00690_FID123364_mc_ep3d_turbine_GRa__sens1.mat'), 'sens');
end
sens100 = sens;


localpath = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/';
webpath = 'https://zenodo.org/record/4743418/files/';

GRa_name = 'meas_MID00690_FID123364_mc_ep3d_turbine_GRa_slice';
GRb_name = 'meas_MID00692_FID123366_mc_ep3d_turbine_GRb_slice';
U8a_name = 'meas_MID00674_FID123348_mc_ep3d_turbine_UNI_8a_slice';
U8b_name = 'meas_MID00676_FID123350_mc_ep3d_turbine_UNI_8b_slice';
U10a_name = 'meas_MID00678_FID123352_mc_ep3d_turbine_UNI_10a_slice';
U10b_name = 'meas_MID00680_FID123354_mc_ep3d_turbine_UNI_10b_slice';
U46a_name = 'meas_MID00682_FID123356_mc_ep3d_turbine_UNI_46a_slice';
U46b_name = 'meas_MID00684_FID123358_mc_ep3d_turbine_UNI_46b_slice';
U55a_name = 'meas_MID00686_FID123360_mc_ep3d_turbine_UNI_55a_slice';
U55b_name = 'meas_MID00688_FID123362_mc_ep3d_turbine_UNI_55b_slice';
S10_46a_name = 'meas_MID00698_FID123372_mc_ep3d_turbine_SILVER_10_46a_slice';
S10_46b_name = 'meas_MID00700_FID123374_mc_ep3d_turbine_SILVER_10_46b_slice';
S8_55a_name = 'meas_MID00694_FID123368_mc_ep3d_turbine_SILVER_8_55a_slice';
S8_55b_name = 'meas_MID00696_FID123370_mc_ep3d_turbine_SILVER_8_55b_slice';

% download data
for slice = 1:16
    if ~exist([localpath GRa_name num2str(slice) '.mat'], 'file')
        websave([localpath GRa_name num2str(slice) '.mat'], [webpath GRa_name num2str(slice) '.mat']);
    end
    if ~exist([localpath GRb_name num2str(slice) '.mat'], 'file')
        websave([localpath GRb_name num2str(slice) '.mat'], [webpath GRb_name num2str(slice) '.mat']);
    end
    if ~exist([localpath U8a_name num2str(slice) '.mat'], 'file')
        websave([localpath U8a_name num2str(slice) '.mat'], [webpath U8a_name num2str(slice) '.mat']);
    end
    if ~exist([localpath U8b_name num2str(slice) '.mat'], 'file')
        websave([localpath U8b_name num2str(slice) '.mat'], [webpath U8b_name num2str(slice) '.mat']);
    end
    if ~exist([localpath U10a_name num2str(slice) '.mat'], 'file')
        websave([localpath U10a_name num2str(slice) '.mat'], [webpath U10a_name num2str(slice) '.mat']);
    end
    if ~exist([localpath U10b_name num2str(slice) '.mat'], 'file')
        websave([localpath U10b_name num2str(slice) '.mat'], [webpath U10b_name num2str(slice) '.mat']);
    end
    if ~exist([localpath U46a_name num2str(slice) '.mat'], 'file')
        websave([localpath U46a_name num2str(slice) '.mat'], [webpath U46a_name num2str(slice) '.mat']);
    end
    if ~exist([localpath U46b_name num2str(slice) '.mat'], 'file')
        websave([localpath U46b_name num2str(slice) '.mat'], [webpath U46b_name num2str(slice) '.mat']);
    end
    if ~exist([localpath U55a_name num2str(slice) '.mat'], 'file')
        websave([localpath U55a_name num2str(slice) '.mat'], [webpath U55a_name num2str(slice) '.mat']);
    end
    if ~exist([localpath U55b_name num2str(slice) '.mat'], 'file')
        websave([localpath U55b_name num2str(slice) '.mat'], [webpath U55b_name num2str(slice) '.mat']);
    end
    if ~exist([localpath S8_55a_name num2str(slice) '.mat'], 'file')
        websave([localpath S8_55a_name num2str(slice) '.mat'], [webpath S8_55a_name num2str(slice) '.mat']);
    end
    if ~exist([localpath S8_55b_name num2str(slice) '.mat'], 'file')
        websave([localpath S8_55b_name num2str(slice) '.mat'], [webpath S8_55b_name num2str(slice) '.mat']);
    end
    if ~exist([localpath S10_46a_name num2str(slice) '.mat'], 'file')
        websave([localpath S10_46a_name num2str(slice) '.mat'], [webpath S10_46a_name num2str(slice) '.mat']);
    end
    if ~exist([localpath S10_46b_name num2str(slice) '.mat'], 'file')
        websave([localpath S10_46b_name num2str(slice) '.mat'], [webpath S10_46b_name num2str(slice) '.mat']);
    end
end

% perform reconstruction (skips if files already exist)
recon_invivo(1:16,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjA/')



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
% TURBINE acquisition in vivo - tSNR resting state linear recon subj B
clear 
close all
%%%%%%%% SETUP files %%%%%%%% 
% insert option to pre-processed data download here for when twix files are not available



%%%%%%%%%%%%%%%%%%%%%%%%
% Subj B
% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20210316_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')% 
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI_slice',0.3, 6016,1:16);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI_slice',1, 6016,1:16);
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20210316_TURBINE/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));

covfile_name = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/covariance_rs_fMRI.mat';
sens_lowres_name = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI__sens0.3.mat';
sens_highres_name = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI__sens1.mat';

if exist(covfile_name, 'file')
    load(covfile_name, 'ncov');
else
    load(websave(covfile_name, 'https://zenodo.org/record/4743764/files/covariance.mat'), 'ncov');
end
    
if exist(sens_lowres_name, 'file')
    load(sens_lowres_name,'sens');
else
    load(websave(sens_lowres_name, 'https://zenodo.org/record/4743764/files/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI__sens0.3.mat'), 'sens');
end
sens030 = sens;

if exist(sens_highres_name, 'file')
    load(sens_highres_name,'sens');
else
    load(websave(sens_highres_name, 'https://zenodo.org/record/4743764/files/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI__sens1.mat'), 'sens');
end
sens100 = sens;


localpath = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/';
webpath = 'https://zenodo.org/record/4743764/files/';

GR_name = 'meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI_slice';
S_name = 'meas_MID00986_FID123660_mc_ep3d_turbine_SILVER_rs_fMRI_slice';

% download data
for slice = 3:13
    if ~exist([localpath GR_name num2str(slice) '.mat'], 'file')
        websave([localpath GR_name num2str(slice) '.mat'], [webpath GR_name num2str(slice) '.mat']);
    end
    if ~exist([localpath S_name num2str(slice) '.mat'], 'file')
        websave([localpath S_name num2str(slice) '.mat'], [webpath S_name num2str(slice) '.mat']);
    end
end

% reconstruct
recon_invivo_fmri(3:13,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')

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
% TURBINE acquisition in vivo - tSNR resting state nonlinear subj B
clear 
close all
%%%%%%%% SETUP files %%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%
% Subj B
% read_in_and_save_TURBINE('/Users/schauman/Documents/FMRIB Employment/Scanning/20210316_TURBINE/', 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')% 
% sens030 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI_slice',0.3, 6016,1:16);
% sens100 = generate_coil_sensitivity_maps_TURBINE('experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI_slice',1, 6016,1:16);
% map = mapVBVD('/Users/schauman/Documents/FMRIB Employment/Scanning/20210316_TURBINE/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI.dat','removeOS',false);
% ncov = cov(reshape(permute(squeeze(map{1}.noise()),[1,3,4,2]),[],32));

covfile_name = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/covariance_rs_fMRI.mat';
sens_lowres_name = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI__sens0.3.mat';
sens_highres_name = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI__sens1.mat';

if exist(covfile_name, 'file')
    load(covfile_name, 'ncov');
else
    load(websave(covfile_name, 'https://zenodo.org/record/4743764/files/covariance.mat'), 'ncov');
end
    
if exist(sens_lowres_name, 'file')
    load(sens_lowres_name,'sens');
else
    load(websave(sens_lowres_name, 'https://zenodo.org/record/4743764/files/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI__sens0.3.mat'), 'sens');
end
sens030 = sens;

if exist(sens_highres_name, 'file')
    load(sens_highres_name,'sens');
else
    load(websave(sens_highres_name, 'https://zenodo.org/record/4743764/files/meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI__sens1.mat'), 'sens');
end
sens100 = sens;


localpath = 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/';
webpath = 'https://zenodo.org/record/4743764/files/';

GR_name = 'meas_MID00988_FID123662_mc_ep3d_turbine_GR_rs_fMRI_slice';
S_name = 'meas_MID00986_FID123660_mc_ep3d_turbine_SILVER_rs_fMRI_slice';

% download data
for slice = 3:13
    if ~exist([localpath GR_name num2str(slice) '.mat'], 'file')
        websave([localpath GR_name num2str(slice) '.mat'], [webpath GR_name num2str(slice) '.mat']);
    end
    if ~exist([localpath S_name num2str(slice) '.mat'], 'file')
        websave([localpath S_name num2str(slice) '.mat'], [webpath S_name num2str(slice) '.mat']);
    end
end

% reconstruct
recon_invivo_fmri_wavelet(3:13,sens030, sens100, 'experiments/data/experiment_inputs/TURBINE_data/in-vivo/subjB/')


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   How much does the trajectory (SILVER, golden ratio, or uniform) affect
%   reconstructed image quality, and does it follow the prediction from
%   earlier efficiency calculations. In-vivo test based on ASL angiography
%   acquired on a 3T Siemens Verio. 
%                                              
%   Sophie Schauman 2020                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% 1. Load data, and choose subject to consider
kdata_GR_filename = 'experiments/data/experiment_inputs/ASL_SILVER/example_kdata_GR_68_153_306.mat';
kdata_S_filename = 'experiments/data/experiment_inputs/ASL_SILVER/example_kdata_SILVER_68_153_306.mat';
kdata_U_filename = 'experiments/data/experiment_inputs/ASL_SILVER/example_kdata_UNIFORM_68_153_306.mat';

sens_filename = 'experiments/data/experiment_inputs/ASL_SILVER/example_sensitivities_68_153_306.mat';

kdata_GR_file = matfile(kdata_GR_filename);
kdata_S_file = matfile(kdata_S_filename);
kdata_U_file = matfile(kdata_U_filename);

sens_file = matfile(sens_filename);


kdata_Uniform = kdata_U_file.kdata_Uniform(:,Subj);
kdata_GR = kdata_GR_file.kdata_GR(:,Subj);
kdata_SILVER = kdata_S_file.kdata_SILVER(:,Subj);

sens_UNIFORM = sens_file.sens_UNIFORM(1,Subj);
sens_GR = sens_file.sens_GR(1,Subj);
sens_SILVER = sens_file.sens_SILVER(1,Subj);


%% 2. Setup useful variables and do SILVER optimiaztion

load('experiments/data/experiment_inputs/ASL_SILVER/example_params_68_153_306.mat')

savename = ['experiments/precalculated_electrostatic_potential/silver_' strrep(num2str(S),' ', '_') '.mat'];
if ~exist(savename, 'file')
    ratio = SILVER_2D(S,'electrostatic_potential',savename) ;
else
    load(savename,'ratio')
end

S_ratio = ratio;

%% 3. Forward transform
savename = ['experiments/data/experiment_results/biorxiv2020/biorxiv2020_in_vivo/subj' num2str(Subj) '/example7_invivo_subj' num2str(Subj) '.mat'];
if ~exist(savename,'file')
    for n = 1:length(S)
        for method = {'UNIFORM', 'GR', 'SILVER'}
            switch method{:}
                case 'UNIFORM'
                    ratio = 1/S(n);
                case 'GR'
                    ratio = gr2D;
                case 'SILVER'
                    ratio = S_ratio;
            end
% OLD METHOD - WRONG!
%             Phi = [];
%             for frame = 1:NFrames
%                 for repeat = 0:NRepeats-1
%                     Phi = cat(1, Phi, mod( ((frame-1)*NSpokes+repeat:NRepeats:frame*NSpokes-1)' * ratio * pi, 2*pi ));
%                 end
%             end
            
            Phi = mod([0:(NSpokes*NFrames-1)] * ratio * pi, 2*pi );
            
            switch method{:}
                case 'UNIFORM'
                    kspace_UNIFORM = reshape(gen_radial_traj(Phi, NSamps, []),[], size(kdata_Uniform{n},2),2);
                case 'GR'
                    kspace_GR = reshape(gen_radial_traj(Phi, NSamps, []),[], size(kdata_Uniform{n},2),2);
                case 'SILVER'
                    kspace_SILVER = reshape(gen_radial_traj(Phi, NSamps, []),[], size(kdata_Uniform{n},2),2);
            end
        end
        E_UNIFORM{n} = xfm_NUFFT([Mat_size,Mat_size,1,size(kdata_Uniform{n},2)],sens_UNIFORM{:},[],...
            kspace_UNIFORM);
        E_GR{n} = xfm_NUFFT([Mat_size,Mat_size,1,size(kdata_Uniform{n},2)],sens_GR{:},[],...
            kspace_GR);
        E_SILVER{n} = xfm_NUFFT([Mat_size,Mat_size,1,size(kdata_Uniform{n},2)],sens_SILVER{:},[],...
            kspace_SILVER);
        clear kspace_UNIFORM
        clear kspace_GR
        clear kspace_SILVER
        clear Phi
    end
        clear sens_UNIFORM
        clear sens_GR
        clear sens_SILVER
else
    warning('using previously reconstructed images')
    load(savename)
end  

%% 3. Reconstruct
if ~exist(savename,'file')
    for n = 1:length(S)
        % linear recon
        kdata_Uniform_presub = kdata_Uniform{n}(:,:,2,:)-kdata_Uniform{n}(:,:,1,:);
        recon_l_Uniform{n} = fista(E_UNIFORM{n}, kdata_Uniform_presub, 1, 0.000000, ...
        [Mat_size, Mat_size, 1,size(kdata_Uniform{n},2)], 100, 0.5);
    
        kdata_GR_presub = kdata_GR{n}(:,:,2,:)-kdata_GR{n}(:,:,1,:);
        recon_l_GR{n} = fista(E_GR{n}, kdata_GR_presub, 1, 0.000000, ...
        [Mat_size, Mat_size, 1,size(kdata_Uniform{n},2)], 100, 0.5);
    
        kdata_SILVER_presub = kdata_SILVER{n}(:,:,2,:)-kdata_SILVER{n}(:,:,1,:);
        recon_l_SILVER{n} = fista(E_SILVER{n}, kdata_SILVER_presub, 1, 0.000000, ...
        [Mat_size, Mat_size, 1,size(kdata_Uniform{n},2)], 100, 0.5);
    
    
        % non-linear recon
        kdata_Uniform_presub = kdata_Uniform{n}(:,:,2,:)-kdata_Uniform{n}(:,:,1,:);
        recon_nl_Uniform{n} = fista(E_UNIFORM{n}, kdata_Uniform_presub, 1, 0.00001, ...
        [Mat_size, Mat_size, 1,size(kdata_Uniform{n},2)], 100, 0.5);
    
        kdata_GR_presub = kdata_GR{n}(:,:,2,:)-kdata_GR{n}(:,:,1,:);
        recon_nl_GR{n} = fista(E_GR{n}, kdata_GR_presub, 1, 0.00001, ...
        [Mat_size, Mat_size, 1,size(kdata_Uniform{n},2)], 100, 0.5);
    
        kdata_SILVER_presub = kdata_SILVER{n}(:,:,2,:)-kdata_SILVER{n}(:,:,1,:);
        recon_nl_SILVER{n} = fista(E_SILVER{n}, kdata_SILVER_presub, 1, 0.00001, ...
        [Mat_size, Mat_size, 1,size(kdata_Uniform{n},2)], 100, 0.5);
    end
    save(savename, 'recon_l_SILVER', 'recon_l_GR', 'recon_l_Uniform', 'recon_nl_SILVER', 'recon_nl_GR', 'recon_nl_Uniform')
end   
 %% 4. Compare reconstructions visually
figure(1)
for n = 1:length(S)
    subplot(length(S),3,n*3-2)
    recon = flipud(abs(mean(recon_l_Uniform{n}(:,:,:,:),4)));
    imagesc(recon)
    cl = caxis;
    caxis([0,0.5*max(cl)])
    axis image
    axis off
    if n ==1
        text(Mat_size/2,-10,'UNIFORM', 'fontsize', 20,'HorizontalAlignment','center')
    end
    text(-2,Mat_size/2, [num2str(S(n)) ' spokes'], 'fontsize', 20,'HorizontalAlignment','right')
    
    
    subplot(length(S),3,n*3-1)
    recon = flipud(abs(mean(recon_l_GR{n}(:,:,:,:),4)));
    imagesc(recon)
    caxis([0,0.5*max(cl)])
    axis image
    axis off
    if n ==1
        text(Mat_size/2,-10,'GOLDEN RATIO', 'fontsize', 20,'HorizontalAlignment','center')
    end
    
    subplot(length(S),3,n*3)
    recon = flipud(abs(mean(recon_l_SILVER{n}(:,:,:,:),4)));
    imagesc(recon)
    caxis([0,0.5*max(cl)])
    axis image
    axis off
    if n ==1
        text(Mat_size/2,-10,'SILVER', 'fontsize', 20,'HorizontalAlignment','center')
    end
    
end
colormap('gray')
set(gcf,'Position', [263 1 827 797])
savefig(['experiments/data/experiment_results/biorxiv2020/subj' num2str(Subj) 'Figure5b.fig'])
saveas(gcf,['experiments/data/experiment_results/biorxiv2020/subj' num2str(Subj) 'Figure5b.tiff'])


figure(2)
for n = 1:length(S)
    subplot(length(S),3,n*3-2)
    recon = flipud(abs(mean(recon_nl_Uniform{n}(:,:,:,:),4)));
    imagesc(recon)
    cl = caxis;
    caxis([0,0.5*max(cl)])
    axis image
    axis off
    if n ==1
        text(Mat_size/2,-10,'UNIFORM', 'fontsize', 20,'HorizontalAlignment','center')
    end
    text(-2,Mat_size/2, [num2str(S(n)) ' spokes'], 'fontsize', 20,'HorizontalAlignment','right')
    
    
    subplot(length(S),3,n*3-1)
    recon = flipud(abs(mean(recon_nl_GR{n}(:,:,:,:),4)));
    imagesc(recon)
    caxis([0,0.5*max(cl)])
    axis image
    axis off
    if n ==1
        text(Mat_size/2,-10,'GOLDEN RATIO', 'fontsize', 20,'HorizontalAlignment','center')
    end
    
    subplot(length(S),3,n*3)
    recon = flipud(abs(mean(recon_nl_SILVER{n}(:,:,:,:),4)));
    imagesc(recon)
    caxis([0,0.5*max(cl)])
    axis image
    axis off
    if n ==1
        text(Mat_size/2,-10,'SILVER', 'fontsize', 20,'HorizontalAlignment','center')
    end
    
end
colormap('gray')
set(gcf,'Position', [263 1 827 797])
savefig(['experiments/data/experiment_results/biorxiv2020/subj' num2str(Subj) '_nonlinear_Figure5b.fig'])
saveas(gcf,['experiments/data/experiment_results/biorxiv2020/subj' num2str(Subj) '_nonlinear_Figure5b.tiff'])

%% 5. Quantify reconstruction quality
% SNR measurement
load('experiments/data/experiment_inputs/ASL_SILVER/example_masks_68_153_306.mat', 'mask_signal', 'mask_noise')

GT = mean(recon_l_Uniform{S == max(S)},4);

for n = 1:length(S)
    S_uniform(n) = mean(abs(recon_l_Uniform{n}(repmat(mask_signal{Subj},[1,1,1,size(recon_l_Uniform{n},4)]))));
    N_uniform(n) = std(abs(recon_l_Uniform{n}(repmat(mask_noise{Subj},[1,1,1,size(recon_l_Uniform{n},4)]))));
    SNR_uniform(n) = S_uniform(n)/N_uniform(n);
    SSIM_uniform(n) = ssim(squeeze(abs(GT./max(GT(:)))),squeeze(abs(mean(recon_nl_Uniform{n},4)./max(mean(recon_nl_Uniform{n},4),[],'all'))));

    
    S_GR(n) = mean(abs(recon_l_GR{n}(repmat(mask_signal{Subj},[1,1,1,size(recon_l_GR{n},4)]))));
    N_GR(n) = std(abs(recon_l_GR{n}(repmat(mask_noise{Subj},[1,1,1,size(recon_l_GR{n},4)]))));
    SNR_GR(n) = S_GR(n)/N_GR(n);
    SSIM_GR(n) = ssim(squeeze(abs(GT./max(GT(:)))),squeeze(abs(mean(recon_nl_GR{n},4)./max(mean(recon_nl_GR{n},4),[],'all'))));

    
    S_SILVER(n) = mean(abs(recon_l_SILVER{n}(repmat(mask_signal{Subj},[1,1,1,size(recon_l_SILVER{n},4)]))));
    N_SILVER(n) = std(abs(recon_l_SILVER{n}(repmat(mask_noise{Subj},[1,1,1,size(recon_l_SILVER{n},4)]))));
    SNR_SILVER(n) = S_SILVER(n)/N_SILVER(n);
    SSIM_SILVER(n) = ssim(squeeze(abs(GT./max(GT(:)))),squeeze(abs(mean(recon_nl_SILVER{n},4)./max(mean(recon_nl_SILVER{n},4),[],'all'))));

    
    
    S_nl_uniform(n) = mean(abs(recon_nl_Uniform{n}(repmat(mask_signal{Subj},[1,1,1,size(recon_nl_Uniform{n},4)]))));
    N_nl_uniform(n) = std(abs(recon_nl_Uniform{n}(repmat(mask_noise{Subj},[1,1,1,size(recon_nl_Uniform{n},4)]))));
    SNR_nl_uniform(n) = S_nl_uniform(n)/N_nl_uniform(n);
    SSIM_l_uniform(n) = ssim(squeeze(abs(GT./max(GT(:)))),squeeze(abs(mean(recon_l_Uniform{n},4)./max(mean(recon_l_Uniform{n},4),[],'all'))));

    
    S_nl_GR(n) = mean(abs(recon_nl_GR{n}(repmat(mask_signal{Subj},[1,1,1,size(recon_nl_GR{n},4)]))));
    N_nl_GR(n) = std(abs(recon_nl_GR{n}(repmat(mask_noise{Subj},[1,1,1,size(recon_nl_GR{n},4)]))));
    SNR_nl_GR(n) = S_nl_GR(n)/N_nl_GR(n);
    SSIM_l_GR(n) = ssim(squeeze(abs(GT./max(GT(:)))),squeeze(abs(mean(recon_l_GR{n},4)./max(mean(recon_l_GR{n},4),[],'all'))));

    
    S_nl_SILVER(n) = mean(abs(recon_nl_SILVER{n}(repmat(mask_signal{Subj},[1,1,1,size(recon_nl_SILVER{n},4)]))));
    N_nl_SILVER(n) = std(abs(recon_nl_SILVER{n}(repmat(mask_noise{Subj},[1,1,1,size(recon_nl_SILVER{n},4)]))));
    SNR_nl_SILVER(n) = S_nl_SILVER(n)/N_nl_SILVER(n);
    SSIM_l_SILVER(n) = ssim(squeeze(abs(GT./max(GT(:)))),squeeze(abs(mean(recon_l_SILVER{n},4)./max(mean(recon_l_SILVER{n},4),[],'all'))));

    
end

figure(3)
subplot(2,2,1)
b = bar(categorical(S), [SNR_uniform; SNR_GR; SNR_SILVER]','BarWidth',1);
b(1).FaceColor = [0,0.5,1];
b(2).FaceColor = [1,0.5,0];
b(3).FaceColor = [0.5,0.5,0.5];
legend('Uniform','GR', 'SILVER', 'Location', 'eastoutside')
title(['Subj ' num2str(Subj) ' - SNR (linear recon)'])
set(gca,'FontSize', 20)
set(gca, 'LineWidth', 2)
grid on
ylabel('SNR')

subplot(2,2,2)
b = bar(categorical(S), [SSIM_uniform; SSIM_GR; SSIM_SILVER]','BarWidth',1);
b(1).FaceColor = [0,0.5,1];
b(2).FaceColor = [1,0.5,0];
b(3).FaceColor = [0.5,0.5,0.5];
legend('Uniform','GR', 'SILVER', 'Location', 'eastoutside')
title(['Subj ' num2str(Subj) ' - SSIM (nonlinear recon)'])
set(gca,'FontSize', 20)
set(gca, 'LineWidth', 2)
grid on
ylabel('SSIM')

subplot(2,2,3)
b = bar(categorical(S), [SNR_nl_uniform; SNR_nl_GR; SNR_nl_SILVER]','BarWidth',1);
b(1).FaceColor = [0,0.5,1];
b(2).FaceColor = [1,0.5,0];
b(3).FaceColor = [0.5,0.5,0.5];
legend('Uniform','GR', 'SILVER', 'Location', 'eastoutside')
title(['Subj ' num2str(Subj) ' - SNR (non-linear recon)'])
set(gca,'FontSize', 20)
set(gca, 'LineWidth', 2)
grid on
ylabel('SNR')

subplot(2,2,4)
b = bar(categorical(S), [SSIM_l_uniform; SSIM_l_GR; SSIM_l_SILVER]','BarWidth',1);
b(1).FaceColor = [0,0.5,1];
b(2).FaceColor = [1,0.5,0];
b(3).FaceColor = [0.5,0.5,0.5];
legend('Uniform','GR', 'SILVER', 'Location', 'eastoutside')
title(['Subj ' num2str(Subj) ' - SSIM (linear recon)'])
set(gca,'FontSize', 20)
set(gca, 'LineWidth', 2)
grid on
ylabel('SSIM')


set(gcf,'Position',[124 359 876 439])

savefig(['experiments/data/experiment_results/biorxiv2020/subj' num2str(Subj) 'Figure5d.fig'])
saveas(gcf,['experiments/data/experiment_results/biorxiv2020/subj' num2str(Subj) 'Figure5.tiff'])


%%
figure(4)
imshow(flipud(abs(mean(recon_l_Uniform{3}(:,:,:,:),4))),[0, 2e-4])
hold on
h = imagesc(cat(3,flipud(mask_signal{Subj}),zeros(192), zeros(192)));
h.AlphaData = 0.2;

hh = imagesc(cat(3,zeros(192), zeros(192), flipud(mask_noise{Subj})));
hh.AlphaData = 0.2;

title(['Subj ' num2str(Subj)], 'Fontsize', 20)
set(gcf, 'Position', [440 84 484 455])

savefig(['experiments/data/experiment_results/biorxiv2020/subj' num2str(Subj) 'FigureS1.fig'])
saveas(gcf,['experiments/data/experiment_results/biorxiv2020/subj' num2str(Subj) 'FigureS1.tiff'])



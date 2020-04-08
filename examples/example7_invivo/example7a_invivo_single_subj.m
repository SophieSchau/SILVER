%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   How much does the trajectory (SILVER, golden ratio, or uniform) affect
%   reconstructed image quality, and does it follow the prediction from
%   earlier efficiency calculations. In-vivo test based on ASL angiography
%   acquired on a 3T Siemens Verio. 
%                                              
%   Sophie Schauman 2020                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% 1. Load data, and choose subject to consider
load('examples/example7_invivo/example_kdata_68_153_306.mat');
load('examples/example7_invivo/example_sensitivities_68_153_306.mat');

Subj = 1;

kdata_Uniform = kdata_Uniform(:,Subj);
kdata_GR = kdata_GR(:,Subj);
kdata_SILVER = kdata_SILVER(:,Subj);

sens_UNIFORM = sens_UNIFORM(Subj);
sens_GR = sens_GR(Subj);
sens_SILVER = sens_SILVER(Subj);


%% 2. Setup useful variables and do SILVER optimiaztion
Mat_size = SILVER_twix.hdr.Config.BaseResolution;
NFrames = SILVER_twix.hdr.Meas.NPhs;
NRepeats = SILVER_twix.hdr.Config.NLin/SILVER_twix.hdr.Config.NSeg;
NSpokes = SILVER_twix.hdr.Config.NLin;
NSamps = SILVER_twix.hdr.Config.NColMeas;
NSubj = size(kdata_SILVER,2);
E_UNIFORM = cell(size(kdata_Uniform));
E_GR = cell(size(kdata_GR));
E_SILVER = cell(size(kdata_SILVER));

clear SILVER_twix;

savename = ['examples/precalculated/silver_' strrep(num2str(S),' ', '_') '.mat'];
if ~exist(savename, 'file')
    ratio = SILVER_2D(S,'electrostatic_potential',savename) ;
else
    load(savename,'ratio')
end

S_ratio = ratio;

%% 3. Forward transform
savename = ['examples/example7_invivo/example7_invivo_subj' num2str(Subj) '.mat'];
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
            Phi = [];
            for frame = 1:NFrames
                for repeat = 0:NRepeats-1
                    Phi = cat(1, Phi, mod( ((frame-1)*NSpokes+repeat:NRepeats:frame*NSpokes-1)' * ratio * pi, 2*pi ));
                end
            end

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
        kdata_Uniform_presub = kdata_Uniform{n}(:,:,2,:)-kdata_Uniform{n}(:,:,1,:);
        recon_Uniform{n} = fista(E_UNIFORM{n}, kdata_Uniform_presub, 1, 0.0000001, ...
        [Mat_size, Mat_size, 1,size(kdata_Uniform{n},2)], 10, 0.5);
    
        kdata_GR_presub = kdata_GR{n}(:,:,2,:)-kdata_GR{n}(:,:,1,:);
        recon_GR{n} = fista(E_GR{n}, kdata_GR_presub, 1, 0.0000001, ...
        [Mat_size, Mat_size, 1,size(kdata_Uniform{n},2)], 10, 0.5);
    
        kdata_SILVER_presub = kdata_SILVER{n}(:,:,2,:)-kdata_SILVER{n}(:,:,1,:);
        recon_SILVER{n} = fista(E_SILVER{n}, kdata_SILVER_presub, 1, 0.0000001, ...
        [Mat_size, Mat_size, 1,size(kdata_Uniform{n},2)], 10, 0.5);
    end
    save(savename, 'recon_SILVER', 'recon_GR', 'recon_Uniform')
end   
 %% 4. Compare reconstructions visually
figure(1)
for n = 1:length(S)
    subplot(length(S),3,n*3-2)
    recon = flipud(abs(mean(recon_Uniform{n}(:,:,:,:),4)));
    imagesc(recon)
    axis image
    axis off
    if n ==1
        text(Mat_size/2,-10,'UNIFORM', 'fontsize', 20,'HorizontalAlignment','center')
    end
    text(-2,Mat_size/2, [num2str(S(n)) ' spokes'], 'fontsize', 20,'HorizontalAlignment','right')
    
    
    subplot(length(S),3,n*3-1)
    recon = flipud(abs(mean(recon_GR{n}(:,:,:,:),4)));
    imagesc(recon)
    axis image
    axis off
    if n ==1
        text(Mat_size/2,-10,'GOLDEN RATIO', 'fontsize', 20,'HorizontalAlignment','center')
    end
    
    subplot(length(S),3,n*3)
    recon = flipud(abs(mean(recon_SILVER{n}(:,:,:,:),4)));
    imagesc(recon)
    axis image
    axis off
    if n ==1
        text(Mat_size/2,-10,'SILVER', 'fontsize', 20,'HorizontalAlignment','center')
    end
    
end

set(gcf,'Position', [263 1 971 797])
savefig(['examples/example7_invivo/example7_invivo_guali_subj_' num2str(Subj) '_result.fig'])
saveas(gcf,['examples/example7_invivo/example7_invivo_guali_subj_' num2str(Subj) '_result.tiff'])

%% 5. Quantify reconstruction quality
% SNR measurement
load('examples/example7_invivo/example_masks_68_153_306.mat', 'mask_signal', 'mask_noise')
for n = 1:length(S)
    S_uniform(n) = mean(abs(recon_Uniform{n}(repmat(mask_signal{Subj},[1,1,1,size(recon_Uniform{n},4)]))));
    N_uniform(n) = std(abs(recon_Uniform{n}(repmat(mask_noise{Subj},[1,1,1,size(recon_Uniform{n},4)]))));
    SNR_uniform(n) = S_uniform(n)/N_uniform(n);
    
    S_GR(n) = mean(abs(recon_GR{n}(repmat(mask_signal{Subj},[1,1,1,size(recon_GR{n},4)]))));
    N_GR(n) = std(abs(recon_GR{n}(repmat(mask_noise{Subj},[1,1,1,size(recon_GR{n},4)]))));
    SNR_GR(n) = S_GR(n)/N_GR(n);
    
    S_SILVER(n) = mean(abs(recon_SILVER{n}(repmat(mask_signal{Subj},[1,1,1,size(recon_SILVER{n},4)]))));
    N_SILVER(n) = std(abs(recon_SILVER{n}(repmat(mask_noise{Subj},[1,1,1,size(recon_SILVER{n},4)]))));
    SNR_SILVER(n) = S_SILVER(n)/N_SILVER(n);

    
end

figure(2)
bar(categorical(S), [SNR_uniform; SNR_GR; SNR_SILVER]','BarWidth',1)
legend('Uniform','GR', 'SILVER', 'Location', 'northwest')
title(['Subj ' num2str(Subj)])
set(gca,'FontSize', 20)
ylabel('SNR')

savefig(['examples/example7_invivo/example7_invivo_SNR_subj_' num2str(Subj) '_result.fig'])
saveas(gcf,['examples/example7_invivo/example7_invivo_SNR_subj_' num2str(Subj) '_result.tiff'])



% % SNR increase of SILVER compared to GR?
% for n = 1:length(S)
%     SNR_prcnt_incr(n) = (SNR_SILVER(n)/SNR_GR(n)-1)*100;
%     theorethical_eff_incr(n) = (efficiency_2D([0:S(n)-1]*ratio*pi, 'Electrostatic_potential')/...
%         efficiency_2D([0:S(n)-1]*gr2D*pi, 'Electrostatic_potential')-1)*100;
% end
% figure(3)
% bar(categorical(S), [SNR_prcnt_incr; theorethical_eff_incr]','BarWidth',1)
% legend('SNR', 'EP', 'Location', 'northeast')
% ylabel('% increase')
% title(['Subj ' num2str(Subj)])
% set(gca,'FontSize', 20)

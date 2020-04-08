%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   How much does the trajectory (SILVER, golden ratio, or uniform) affect
%   reconstructed image quality, and does it follow the prediction from
%   earlier efficiency calculations. Phantom test.
%                                              
%   Sophie Schauman 2020                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close 

%% 1. Choose  set of window sizes, S, to consider
S = [16,32,64];

%% 2. Do SILVER optimization for that range

savename = ['examples/precalculated/silver_' strrep(num2str(S),' ', '_') '.mat'];
if ~exist(savename, 'file')
    ratio = SILVER_2D(S,'electrostatic_potential',savename) ;
else
    load(savename,'ratio')
end

%% 3. Set up ground truth and forward model
load('examples/example6_phantom_simu/phantom_64x64x1x60','im');
rng(123)

for n = 1:length(S)
    k_Uniform= reshape(gen_radial_traj((0:60*S(n)-1)*1/S(n)*pi, 128, []),[],60,2);
    k_GR= reshape(gen_radial_traj((0:60*S(n)-1)*gr2D*pi, 128, []),[],60,2);
    k_SILVER= reshape(gen_radial_traj((0:60*S(n)-1)*ratio*pi, 128, []),[],60,2);
    
    E_Uniform{n} = xfm_NUFFT([64 64 1 60],ones(64),[],k_Uniform);
    E_GR{n} = xfm_NUFFT([64 64 1 60],ones(64),[],k_GR);
    E_SILVER{n} = xfm_NUFFT([64 64 1 60],ones(64),[],k_SILVER);
    
    noise = (randn(size(k_Uniform,1),size(k_Uniform,2)) + 1i* randn(size(k_Uniform,1),size(k_Uniform,2)))/sqrt(2);
    
    kdata_Uniform{n} = E_Uniform{n}*im +noise;
    kdata_GR{n} = E_GR{n}*im +noise;
    kdata_SILVER{n} = E_SILVER{n}*im +noise;
end

%% 4. Reconstruct
savename = 'examples/example6_phantom_simu/example6_phantom_simu.mat';
if ~exist(savename,'file')
    for n = 1:length(S)
        recon_Uniform{n} = cham_pock_2_xtTGV(E_Uniform{n}'*kdata_Uniform{n}, E_Uniform{n}, 500, 0.5, 0);
        recon_GR{n} = cham_pock_2_xtTGV(E_GR{n}'*kdata_GR{n}, E_GR{n}, 500, 0.5, 0);
        recon_SILVER{n} = cham_pock_2_xtTGV(E_SILVER{n}'*kdata_SILVER{n}, E_SILVER{n}, 500, 0.5, 0);
    end
    save(savename, 'recon_SILVER', 'recon_GR', 'recon_Uniform')
else
    warning('using previously reconstructed images')
    load(savename)
end   
%% 5. Compare reconstructions visually
figure(1)
for n = 1:length(S)
    subplot(length(S),3,n*3-2)
    recon = abs(recon_Uniform{n}.out(:,:,:,10));
    imagesc(recon)
    axis image
    axis off
    if n ==1
        text(32,-4,'UNIFORM', 'fontsize', 20,'HorizontalAlignment','center')
    end
    text(-2,32, [num2str(S(n)) ' spokes'], 'fontsize', 20,'HorizontalAlignment','right')
    
    
    subplot(length(S),3,n*3-1)
    recon = abs(recon_GR{n}.out(:,:,:,10));
    imagesc(recon)
    axis image
    axis off
    if n ==1
        text(32,-4,'GOLDEN RATIO', 'fontsize', 20,'HorizontalAlignment','center')
    end
    
    subplot(length(S),3,n*3)
    recon = abs(recon_SILVER{n}.out(:,:,:,10));
    imagesc(recon)
    axis image
    axis off
    if n ==1
        text(32,-4,'SILVER', 'fontsize', 20,'HorizontalAlignment','center')
    end
    
end

set(gcf,'Position', [440 1 794 797])
savefig('examples/example6_phantom_simu/example6_phantom_simu_guali_result.fig')
saveas(gcf,'examples/example6_phantom_simu/example6_phantom_simu_quali_result.tiff')

%% 6. Quantify reconstruction quality
% SNR and SSIM measurement
signal_mask = logical(im);
noise_mask = ~signal_mask;
for n = 1:length(S)
    S_uniform(n) = mean(abs(recon_Uniform{n}.out(signal_mask)));
    N_uniform(n) = std(recon_Uniform{n}.out(noise_mask));
    SNR_uniform(n) = S_uniform(n)/N_uniform(n);
    SSIM_uniform(n) = ssim(single(squeeze(abs(im./max(im(:))))),squeeze(abs(recon_Uniform{n}.out./max(recon_Uniform{n}.out(:)))));
    
    S_GR(n) = mean(abs(recon_GR{n}.out(signal_mask)));
    N_GR(n) = std(recon_GR{n}.out(noise_mask));
    SNR_GR(n) = S_GR(n)/N_GR(n);
    SSIM_GR(n) = ssim(single(squeeze(abs(im./max(im(:))))),squeeze(abs(recon_GR{n}.out./max(recon_GR{n}.out(:)))));

    
    S_SILVER(n) = mean(abs(recon_SILVER{n}.out(signal_mask)));
    N_SILVER(n) = std(recon_SILVER{n}.out(noise_mask));
    SNR_SILVER(n) = S_SILVER(n)/N_SILVER(n);
    SSIM_SILVER(n) = ssim(single(squeeze(abs(im./max(im(:))))),squeeze(abs(recon_SILVER{n}.out./max(recon_SILVER{n}.out(:)))));

    
end

figure(2)
bar(categorical(S), [SNR_uniform; SNR_GR; SNR_SILVER]')
legend('Uniform','GR', 'SILVER', 'Location', 'northwest')
set(gca,'FontSize', 20)
ylabel('SNR')

savefig('examples/example6_phantom_simu/example6_phantom_simu_SNR_result.fig')
saveas(gcf,'examples/example6_phantom_simu/example6_phantom_simu_SNR_result.tiff')


figure(3)
bar(categorical(S), [SSIM_uniform; SSIM_GR; SSIM_SILVER]')
legend('Uniform','GR', 'SILVER', 'Location', 'northwest')
set(gca,'FontSize', 20)
ylabel('SSIM')

savefig('examples/example6_phantom_simu/example6_phantom_simu_SSIM_result.fig')
saveas(gcf,'examples/example6_phantom_simu/example6_phantom_simu_SSIM_result.tiff')

% SNR and SSIM increase of SILVER compared to GR?
for n = 1:length(S)
    SNR_prcnt_incr(n) = (SNR_SILVER(n)/SNR_GR(n)-1)*100;
    SSIM_prcnt_incr(n) = (SSIM_SILVER(n)/SSIM_GR(n)-1)*100;
    theorethical_eff_incr(n) = (efficiency_2D([0:S(n)-1]*ratio*pi, 'Electrostatic_potential')/...
        efficiency_2D([0:S(n)-1]*gr2D*pi, 'Electrostatic_potential')-1)*100;
end
figure(4)
bar(categorical(S), [SNR_prcnt_incr; SSIM_prcnt_incr; theorethical_eff_incr]')
legend('SNR','SSIM', 'EP', 'Location', 'northeast')
ylabel('% increase')
set(gca,'FontSize', 20)
savefig('examples/example6_phantom_simu/example6_phantom_simu_quality_incr_result.fig')
saveas(gcf,'examples/example6_phantom_simu/example6_phantom_simu_quality_incr_result.tiff')
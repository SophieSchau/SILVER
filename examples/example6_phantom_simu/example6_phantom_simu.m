%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   How much does the trajectory (SILVER, golden ratio, or uniform) affect
%   reconstructed image quality, and does it follow the prediction from
%   earlier efficiency calculations. Phantom test.
%                                              
%   Sophie Schauman 2020                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% 1. Choose  set of window sizes, S, to consider
S = [16,32,48];

%% 2. Do SILVER optimization for that range

savename = ['examples/precalculated/silver_' strrep(num2str(S),' ', '_') '.mat'];
if ~exist(savename, 'file')
    ratio = SILVER_2D(S,'electrostatic_potential',savename) ;
else
    load(savename,'ratio')
end

%% 3. Set up ground truth, forward model, and reconstruct
load('examples/example6_phantom_simu/phantom_64x64x1x60','im');
load('examples/example5_gfactor/test_sensitivity_map.mat','psens');


savename = 'examples/example6_phantom_simu/example6_phantom_simu.mat';

if ~exist(savename,'file')
    
    for n = 1:length(S)
        k_Uniform= reshape(gen_radial_traj((0:60*S(n)-1)*1/S(n)*pi, 128, []),[],60,2);
        k_GR= reshape(gen_radial_traj((0:60*S(n)-1)*gr2D*pi, 128, []),[],60,2);
        k_SILVER= reshape(gen_radial_traj((0:60*S(n)-1)*ratio*pi, 128, []),[],60,2);

        E_Uniform{n} = xfm_NUFFT([64 64 1 60],psens,[],k_Uniform);
        E_GR{n} = xfm_NUFFT([64 64 1 60],psens,[],k_GR);
        E_SILVER{n} = xfm_NUFFT([64 64 1 60],psens,[],k_SILVER);
        for seed = 1:10
            rng(seed)
            noise = (randn(size(k_Uniform,1),size(k_Uniform,2)) + 1i* randn(size(k_Uniform,1),size(k_Uniform,2)))/sqrt(2);
            kdata_Uniform = E_Uniform{n}*im +noise;
            kdata_GR = E_GR{n}*im +noise;
            kdata_SILVER = E_SILVER{n}*im +noise;

            recon_l_Uniform{n,seed} = fista(E_Uniform{n}, kdata_Uniform, 1, 0.000000, ...
            [64, 64, 1,size(kdata_Uniform,2)], 50, 0.5);
            recon_l_GR{n,seed} = fista(E_GR{n}, kdata_GR, 1, 0.0000000, ...
            [64, 64, 1,size(kdata_GR,2)], 50, 0.5);
            recon_l_SILVER{n,seed} = fista(E_SILVER{n}, kdata_SILVER, 1, 0.000000, ...
            [64, 64, 1,size(kdata_SILVER,2)], 50, 0.5);
        end
    end
    save(savename, 'recon_l_SILVER', 'recon_l_GR', 'recon_l_Uniform')

else
    warning('using previously reconstructed images')
    load(savename, 'recon_l_SILVER', 'recon_l_GR', 'recon_l_Uniform')
end

%% 4. Compare reconstructions visually
figure(1)
clim = [min(im(:)),max(im(:))];
for t = 1:60
for n = 1:length(S)
    subplot(length(S),3,n*3-2)
    recons = [recon_l_Uniform{n,1}];
    recons = permute(reshape(recons, 64,64,1,1,60),[1,2,4,5,3]);
    recons_mean = mean(recons,5);
    imagesc(cat(2,abs(recons_mean(:,:,1,t)),abs(im(:,:,1,t)-recons_mean(:,:,1,t))) );
    axis image
    axis off
%     caxis(clim)
    if n ==1
        text(64,-6,'UNIFORM', 'fontsize', 20,'HorizontalAlignment','center')
    end
    text(-2,32, [num2str(S(n)) ' spokes'], 'fontsize', 20,'HorizontalAlignment','right')
    
    
    subplot(length(S),3,n*3-1)
    recons = [recon_l_GR{n,1}];
    recons = permute(reshape(recons, 64,64,1,1,60),[1,2,4,5,3]);
    recons_mean = mean(recons,5);
    imagesc(cat(2,abs(recons_mean(:,:,1,t)),abs(im(:,:,1,t)-recons_mean(:,:,1,t))) );
    axis image
    axis off
%     caxis(clim)
    if n ==1
        text(64,-6,'GOLDEN RATIO', 'fontsize', 20,'HorizontalAlignment','center')
    end
    
    subplot(length(S),3,n*3)
    recons = [recon_l_SILVER{n,1}];
    recons = permute(reshape(recons, 64,64,1,1,60),[1,2,4,5,3]);
    recons_mean = mean(recons,5);
    imagesc(cat(2,abs(recons_mean(:,:,1,t)),abs(im(:,:,1,t)-recons_mean(:,:,1,t))) );
    axis image
    axis off
%     caxis(clim)
    if n ==1
        text(64,-6,'SILVER', 'fontsize', 20,'HorizontalAlignment','center')
    end
    
end
drawnow
set(gcf,'Position', [81 230 1001 490])
makegif('examples/example6_phantom_simu/example6_phantom_simu_guali_result.gif')
end
savefig('examples/example6_phantom_simu/example6_phantom_simu_guali_result.fig')
saveas(gcf,'examples/example6_phantom_simu/example6_phantom_simu_quali_result.tiff')

%% 6. Quantify reconstruction quality
% SNR and SSIM measurement
signal_mask = logical(im);
noise_mask = ~signal_mask;
for n = 1:length(S)
    for seed = 1:10
        S_uniform(n,seed) = mean(abs(recon_l_Uniform{n,seed}(signal_mask)));
        N_uniform(n,seed) = std(recon_l_Uniform{n,seed}(noise_mask));
        SNR_uniform(n,seed) = S_uniform(n,seed)/N_uniform(n,seed);

        S_GR(n,seed) = mean(abs(recon_l_GR{n,seed}(signal_mask)));
        N_GR(n,seed) = std(recon_l_GR{n,seed}(noise_mask));
        SNR_GR(n,seed) = S_GR(n,seed)/N_GR(n,seed);

        S_SILVER(n,seed) = mean(abs(recon_l_SILVER{n,seed}(signal_mask)));
        N_SILVER(n,seed) = std(recon_l_SILVER{n,seed}(noise_mask));
        SNR_SILVER(n,seed) = S_SILVER(n,seed)/N_SILVER(n,seed);
    end
    
end

figure(2)
res = [mean(SNR_uniform,2,'omitnan' )'; mean(SNR_GR,2,'omitnan' )'; mean(SNR_SILVER,2,'omitnan' )']';
err = [std(SNR_uniform,0,2,'omitnan' )'; std(SNR_GR,0,2,'omitnan' )'; std(SNR_SILVER,0,2,'omitnan' )']';

bar(res, 'barwidth', 1)
hold on
legend('Uniform','GR', 'SILVER', 'Location', 'northwest')
set(gca,'FontSize', 20)
ylabel('SNR')
xticklabels(S)

ngroups = length(S);
nbars = 3;
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x,res(:,i), err(:,i),'k', 'linestyle', 'none', 'linewidth', 1);    
end
x = (1:ngroups) - groupwidth/2 + groupwidth / (2*nbars);
yt = get(gca, 'YTick');
axis([xlim    0  ceil(max(yt)*1.15)])
for s = 1:length(S)
    if ttest2(SNR_uniform(s,:),SNR_GR(s,:))
        plot([x(s),x(s)+groupwidth/3], [1 1]*max(res(s,:))+2.2, '-k',  mean([x(s),x(s)+groupwidth/3]), max(res(s,:))+2.3, '*k')
    end
    if ttest2(SNR_uniform(s,:),SNR_SILVER(s,:))
        yt = get(gca, 'YTick');
        axis([xlim    0  ceil(max(yt)*1.15)])
        hold on
        plot([x(s),x(s)+groupwidth*2/3], [1 1]*max(res(s,:))+1.2, '-k',  mean([x(s),x(s)+groupwidth*2/3]), max(res(s,:))+1.3, '*k')
    end
    if ttest2(SNR_SILVER(s,:),SNR_GR(s,:))
        yt = get(gca, 'YTick');
        axis([xlim    0  ceil(max(yt)*1.15)])
        hold on
        plot([x(s)+groupwidth/3,x(s)+groupwidth*2/3], [1 1]*max(res(s,:))+3.2, '-k',  mean([x(s)+groupwidth/3,x(s)+groupwidth*2/3]), max(res(s,:))+3.3, '*k')
    end
end


legend('Uniform','GR', 'SILVER', 'Location', 'northwest')

savefig('examples/example6_phantom_simu/example6_phantom_simu_SNR_result.fig')
saveas(gcf,'examples/example6_phantom_simu/example6_phantom_simu_SNR_result.tiff')
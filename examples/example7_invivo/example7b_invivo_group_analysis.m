%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   How much does the trajectory (SILVER, golden ratio, or uniform) affect
%   reconstructed image quality, and does it follow the prediction from
%   earlier efficiency calculations. In-vivo test based on ASL angiography
%   acquired on a 3T Siemens Verio. Comparing multiple subjects.
%                                              
%   Sophie Schauman 2020                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
close all

%% Group analysis

subjects = [1,2,3];
S = [68,153,306];
load('examples/example7_invivo/example_masks_68_153_306.mat', 'mask_signal', 'mask_noise')

for subj = subjects
    savename = ['examples/example7_invivo/subj' num2str(subj) '/example7_invivo_subj' num2str(subj) '.mat'];
    load(savename, 'recon_l_SILVER', 'recon_l_GR', 'recon_l_Uniform')
    for n = 1:length(S)
        S_uniform(n,subj) = mean(abs(recon_l_Uniform{n}(repmat(mask_signal{subj},[1,1,1,size(recon_l_Uniform{n},4)]))));
        N_uniform(n,subj) = std(abs(recon_l_Uniform{n}(repmat(mask_noise{subj},[1,1,1,size(recon_l_Uniform{n},4)]))));
        SNR_uniform(n,subj) = S_uniform(n,subj)/N_uniform(n,subj);

        S_GR(n,subj) = mean(abs(recon_l_GR{n}(repmat(mask_signal{subj},[1,1,1,size(recon_l_GR{n},4)]))));
        N_GR(n,subj) = std(abs(recon_l_GR{n}(repmat(mask_noise{subj},[1,1,1,size(recon_l_GR{n},4)]))));
        SNR_GR(n,subj) = S_GR(n,subj)/N_GR(n,subj);

        S_SILVER(n,subj) = mean(abs(recon_l_SILVER{n}(repmat(mask_signal{subj},[1,1,1,size(recon_l_SILVER{n},4)]))));
        N_SILVER(n,subj) = std(abs(recon_l_SILVER{n}(repmat(mask_noise{subj},[1,1,1,size(recon_l_SILVER{n},4)]))));

        SNR_SILVER(n,subj) = S_SILVER(n,subj)/N_SILVER(n,subj);
    end
end
SNR_uniform(SNR_uniform==0)=NaN;
SNR_GR(SNR_GR==0)=NaN;
SNR_SILVER(SNR_SILVER==0)=NaN;


res = [mean(SNR_uniform,2,'omitnan' )'; mean(SNR_GR,2,'omitnan' )'; mean(SNR_SILVER,2,'omitnan' )']';
err = [std(SNR_uniform,0,2,'omitnan' )'; std(SNR_GR,0,2,'omitnan' )'; std(SNR_SILVER,0,2,'omitnan' )']';

figure(4)
bar(res,'BarWidth',1)
xticklabels(S)
xlabel('Number of spokes')
set(gca,'FontSize', 20)
title(['Subjects ' num2str(subjects)])
ylabel('SNR')
hold on

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
axis([xlim    0  ceil(max(yt)*1.05)])
for s = 1:length(S)
    if ttest2(SNR_uniform(s,:),SNR_GR(s,:))
        plot([x(s),x(s)+groupwidth/3], [1 1]*max(res(s,:))*1.1, '-k',  mean([x(s),x(s)+groupwidth/3]), max(res(s,:))*1.15, '*k')
    end
    if ttest2(SNR_uniform(s,:),SNR_SILVER(s,:))
        yt = get(gca, 'YTick');
        axis([xlim    0  ceil(max(yt)*1.05)])
        hold on
        plot([x(s),x(s)+groupwidth*2/3], [1 1]*max(res(s,:))*1.1, '-k',  mean([x(s),x(s)+groupwidth*2/3]), max(res(s,:))*1.15, '*k')
    end
    if ttest2(SNR_SILVER(s,:),SNR_GR(s,:))
        yt = get(gca, 'YTick');
        axis([xlim    0  ceil(max(yt)*1.05)])
        hold on
        plot([x(s)+groupwidth/3,x(s)+groupwidth*2/3], [1 1]*max(res(s,:))*1.1, '-k',  mean([x(s)+groupwidth/3,x(s)+groupwidth*2/3]), max(res(s,:))*1.15, '*k')
    end
end
        

legend('Uniform','GR', 'SILVER', 'Location', 'northwest')

savefig(['examples/example7_invivo/group/example7_invivo_SNR_group_result.fig'])
saveas(gcf,['examples/example7_invivo/group/example7_invivo_SNR_group_result.tiff'])


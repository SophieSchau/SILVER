function [] = tSNR_fig(recons_folder,slices_for_calc, slice_for_viz, mask_lowres, mask_highres, savename)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

for sl = slices_for_calc
    GR_file_10 = dir([recons_folder 'Recon_sl' num2str(sl) '_ws10*GR_rs_*']);
    GR_file_46 = dir([recons_folder 'Recon_sl' num2str(sl) '_ws46*GR_rs_*']);
    S_file_10 = dir([recons_folder 'Recon_sl' num2str(sl) '_ws10*SILVER_rs_*']);
    S_file_46 = dir([recons_folder 'Recon_sl' num2str(sl) '_ws46*SILVER_rs_*']);
    
    load([recons_folder GR_file_10.name], 'recon_a')
    GR_10(:,:,sl,:) = abs(recon_a);
    
    load([recons_folder GR_file_46.name], 'recon_a')
    GR_46(:,:,sl,:) = abs(recon_a);
    
    load([recons_folder S_file_10.name], 'recon_a')
    S_10(:,:,sl,:) = abs(recon_a);
    
    load([recons_folder S_file_46.name], 'recon_a')
    S_46(:,:,sl,:) = abs(recon_a);
end
GR_rs_S_10 = mean(GR_10,4);
GR_rs_N_10 = std(GR_10,[],4);
tSNR_GR_10 = GR_rs_S_10./GR_rs_N_10;

GR_rs_S_46 = mean(GR_46,4);
GR_rs_N_46 = std(GR_46,[],4);
tSNR_GR_46 = GR_rs_S_46./GR_rs_N_46;

SILVER_rs_S_10 = mean(S_10,4);
SILVER_rs_N_10 = std(S_10,[],4);
tSNR_SILVER_10 = SILVER_rs_S_10./SILVER_rs_N_10;

SILVER_rs_S_46 = mean(S_46,4);
SILVER_rs_N_46 = std(S_46,[],4);
tSNR_SILVER_46 = SILVER_rs_S_46./SILVER_rs_N_46;

subplot(2,3,1)
imagesc(imrotate(tSNR_GR_10(:,:,slice_for_viz).*mask_lowres(:,:,slice_for_viz),-90))
tmpA = tSNR_GR_10(:,:,slices_for_calc);
text(2,2,['tSNR = ' num2str(mean(tmpA(mask_lowres(:,:,slices_for_calc)),'all'))], 'FontSize',16, 'Color','w')
axis image
axis off
colormap hot
caxis([0,50])
colorbar

text(-5, 15, 'Navigator - 10 spokes/frame', 'Rotation', 90, 'FontSize', 16, 'HorizontalAlignment', 'center');
text(15,-3, 'Golden ratio', 'FontSize',16, 'HorizontalAlignment', 'center'); 
text(-5*2,-3, '(A)', 'FontSize',16, 'HorizontalAlignment', 'center'); 


subplot(2,3,2)
imagesc(imrotate(tSNR_SILVER_10(:,:,slice_for_viz).*mask_lowres(:,:,slice_for_viz),-90))
tmpB = tSNR_SILVER_10(:,:,slices_for_calc);
text(2,2,['tSNR = ' num2str(mean(tmpB(mask_lowres(:,:,slices_for_calc)),'all'))], 'FontSize',16, 'Color', 'w')
axis image
axis off
colormap hot
caxis([0,50])
colorbar

text(15,-3, 'SILVER', 'FontSize',16, 'HorizontalAlignment', 'center'); 



subplot(2,3,4)
imagesc(imrotate(tSNR_GR_46(:,:,slice_for_viz).*mask_highres(:,:,slice_for_viz),-90))
tmpC = tSNR_GR_46(:,:,slices_for_calc);
text(7,7,['tSNR = ' num2str(mean(tmpC(mask_highres(:,:,slices_for_calc)),'all'))], 'FontSize',16, 'Color', 'w')
axis image
axis off
colormap hot
caxis([0,50])
colorbar

text(50,-10, 'Golden ratio', 'FontSize',16, 'HorizontalAlignment', 'center'); 
text(-17, 50, 'Imaging - 46 spokes/frame', 'Rotation', 90, 'FontSize', 16, 'HorizontalAlignment', 'center');
text(-17*2,-10, '(B)', 'FontSize',16, 'HorizontalAlignment', 'center'); 

subplot(2,3,5)
imagesc(imrotate(tSNR_SILVER_46(:,:,slice_for_viz).*mask_highres(:,:,slice_for_viz),-90))
tmpD = tSNR_SILVER_46(:,:,slices_for_calc);
text(7,7,['tSNR = ' num2str(mean(tmpD(mask_highres(:,:,slices_for_calc)),'all'))], 'FontSize',16, 'Color', 'w')
axis image
axis off
colormap hot
caxis([0,50])
colorbar
text(50,-10, 'SILVER', 'FontSize',16, 'HorizontalAlignment', 'center'); 



subplot(2,3,3)
[p, h] = signrank(reshape(tmpA(mask_lowres(:,:,slices_for_calc)),[],1),reshape(tmpB(mask_lowres(:,:,slices_for_calc)),[],1), 'alpha', 0.05/2);

hA = histogram(reshape(tmpA(mask_lowres(:,:,slices_for_calc)),[],1), 'BinWidth',1, 'FaceColor', [1 0.5 0]);
hold on
hB = histogram(reshape(tmpB(mask_lowres(:,:,slices_for_calc)),[],1), 'BinWidth',1, 'FaceColor', [0.5 0.5 0.5]);

if ~h
    [max_a] = max(hA.Values);
    [max_b] = max(hB.Values);
    line([hA.BinEdges(end/2), hB.BinEdges(end/2)], max([max_a,max_b])*[1.001 1.001], 'handlevisibility', 'off', 'color', 'k')
    text(mean([hA.BinEdges(end/2), hB.BinEdges(end/2)]), max([max_a,max_b])*1.01, 'n.s.', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
end

legend('GR', 'SILVER')
title('tSNR - low res.', 'FontSize',16)
axis([0 70 ylim])
xlabel('tSNR')

subplot(2,3,6)
[p, h] = signrank(reshape(tmpC(mask_highres(:,:,slices_for_calc)),[],1),reshape(tmpD(mask_highres(:,:,slices_for_calc)),[],1), 'alpha', 0.05/2);


hC = histogram(reshape(tmpC(mask_highres(:,:,slices_for_calc)),[],1), 'BinWidth',1, 'FaceColor', [1 0.5 0]);
hold on
hD = histogram(reshape(tmpD(mask_highres(:,:,slices_for_calc)),[],1), 'BinWidth',1, 'FaceColor', [0.5 0.5 0.5]);

if ~h
    [max_c] = max(hC.Values);
    [max_d] = max(hD.Values);
    line([hC.BinEdges(end/2), hD.BinEdges(end/2)], max([max_c,max_d])*[1.001 1.001], 'handlevisibility', 'off', 'color', 'k')
    text(mean([hC.BinEdges(end/2), hD.BinEdges(end/2)]), max([max_c,max_d])*1.01, 'n.s.', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
end
legend('GR', 'SILVER')
title('tSNR - high res.', 'FontSize',16)
axis([0 70 ylim])
xlabel('tSNR')

annotation('line',[0 1], [0.5 0.5]);


set(gcf,'Position',[21 81 1242 717])

set(gcf, 'InvertHardCopy', 'off');
set(gcf,'Color', 'w');
savefig(gcf, [savename '.fig'])
saveas(gcf, [savename '.svg'])
end


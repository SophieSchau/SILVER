function [] = rs_recon_fig(recons_folder, slice, titl, savename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

filename_GR_lo = dir([recons_folder 'Recon_sl' num2str(slice) '_ws10*GR_rs_*.mat']);
filename_SILVER_lo = dir([recons_folder 'Recon_sl' num2str(slice) '_ws10*SILVER_rs_*.mat']);
filename_GR_hi = dir([recons_folder 'Recon_sl' num2str(slice) '_ws46*GR_rs_*.mat']);
filename_SILVER_hi = dir([recons_folder 'Recon_sl' num2str(slice) '_ws46*SILVER_rs_*.mat']);


load([recons_folder filename_GR_lo.name], 'recon_a')
GR_lo = recon_a;
clear recon_a;

load([recons_folder filename_SILVER_lo.name], 'recon_a')
SILVER_lo = recon_a;
clear recon_a;

load([recons_folder filename_GR_hi.name], 'recon_a')
GR_hi = recon_a;
clear recon_a;

load([recons_folder filename_SILVER_hi.name], 'recon_a')
SILVER_hi = recon_a;
clear recon_a;


subplot(2,2,1)
imagesc(cat(2,imrotate(mean(abs(GR_lo(:,:,:,:)),4),-90), 10*imrotate(std(abs(GR_lo(:,:,:,:)),[],4),-90)));
axis image
axis off
colormap gray
cl = caxis;
text(15,3, 'GR - mean', 'Color', 'w', 'fontsize', 14, 'HorizontalAlignment', 'center')
text(45,3, 'GR - std x 10', 'Color', 'w', 'fontsize', 14, 'HorizontalAlignment', 'center')
text(-5, 15, 'Navigator - 10 spokes/frame', 'Rotation', 90, 'FontSize', 16, 'HorizontalAlignment', 'center');


subplot(2,2,2)
imagesc(cat(2,imrotate(mean(abs(SILVER_lo(:,:,:,:)),4),-90), 10*imrotate(std(abs(SILVER_lo(:,:,:,:)),[],4),-90)));
axis image
axis off
colormap gray
caxis(cl);
text(15,3, 'SILVER - mean', 'Color', 'w', 'fontsize', 14, 'HorizontalAlignment', 'center')
text(45,3, 'SILVER - std x 10', 'Color', 'w', 'fontsize', 14, 'HorizontalAlignment', 'center')

subplot(2,2,3)
imagesc(cat(2,imrotate(mean(abs(GR_hi(:,:,:,:)),4),-90), 10*imrotate(std(abs(GR_hi(:,:,:,:)),[],4),-90)));
axis image
axis off
colormap gray
cl = caxis;
text(15/0.3,3/0.3, 'GR - mean', 'Color', 'w', 'fontsize', 14, 'HorizontalAlignment', 'center')
text(45/0.3,3/0.3, 'GR - std x 10', 'Color', 'w', 'fontsize', 14, 'HorizontalAlignment', 'center')
text(-17, 50, 'Imaging - 46 spokes/frame', 'Rotation', 90, 'FontSize', 16, 'HorizontalAlignment', 'center');


subplot(2,2,4)
imagesc(cat(2,imrotate(mean(abs(SILVER_hi(:,:,:,:)),4),-90), 10*imrotate(std(abs(SILVER_hi(:,:,:,:)),[],4),-90)));
axis image
axis off
colormap gray
caxis(cl);
text(15/0.3,3/0.3, 'SILVER - mean', 'Color', 'w', 'fontsize', 14, 'HorizontalAlignment', 'center')
text(45/0.3,3/0.3, 'SILVER - std x 10', 'Color', 'w', 'fontsize', 14, 'HorizontalAlignment', 'center')

set(gcf, 'position', [128 42 1263 756])
annotation('textbox', [0.45, 0.9 0.1, 0.1], 'string', titl, 'LineStyle', 'none', 'fontsize', 18, 'HorizontalAlignment', 'center')

set(gcf, 'InvertHardcopy', 'off')
set(gcf, 'color', 'w')
savefig(gcf, [savename '.fig'])
saveas(gcf, [savename '.svg'])
end


clear
close all
for subj = 1:3
    savename = ['examples/example7_invivo/subj' num2str(subj) '/example7_invivo_subj' num2str(subj) 'test.mat'];
    subj_data{subj} = load(savename, 'recon_l_SILVER', 'recon_l_GR', 'recon_l_Uniform');    
end        
S = [68,153,306];
for n = 1:3
    figure
        for t = 1:size(subj_data{1}.recon_l_Uniform{n},4)
            im = cat(1,...
                flipud(abs(cat(2,subj_data{1}.recon_l_Uniform{n}(:,:,:,t),subj_data{1}.recon_l_GR{n}(:,:,:,t), subj_data{1}.recon_l_SILVER{n}(:,:,:,t)))),...
                flipud(abs(cat(2,subj_data{2}.recon_l_Uniform{n}(:,:,:,t),subj_data{2}.recon_l_GR{n}(:,:,:,t), subj_data{2}.recon_l_SILVER{n}(:,:,:,t)))),...
                flipud(abs(cat(2,subj_data{3}.recon_l_Uniform{n}(:,:,:,t),subj_data{3}.recon_l_GR{n}(:,:,:,t), subj_data{3}.recon_l_SILVER{n}(:,:,:,t)))));
           
            
            
            imagesc(im,[0 3e-4])
            text(96,590,'UNIFORM', 'fontsize', 20,'HorizontalAlignment','center')
            text(288,590,'Golden ratio', 'fontsize', 20,'HorizontalAlignment','center')
            text(480,590,'SILVER', 'fontsize', 20,'HorizontalAlignment','center')
            
            text(-10,96,'Subj 1', 'fontsize', 20,'HorizontalAlignment','center', 'rotation', 90)
            text(-10,288,'Subj 2', 'fontsize', 20,'HorizontalAlignment','center','rotation', 90)
            text(-10,480,'Subj 3', 'fontsize', 20,'HorizontalAlignment','center','rotation', 90)
            
            axis image
            axis off
            colormap gray
            title([num2str(S(n)) ' spokes'],'Fontsize', 20)
            set(gcf, 'Position', [20 58 903 740])
            drawnow
            makegif_fast(['examples/example7_invivo/group/example7_invivo_N' num2str(S(n)) '.gif'])
        end
end
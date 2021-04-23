function [outputArg1,outputArg2] = SNR_phantom_quali_fig(sensmask_lowres, sensmask_highres, slices, kdata_folder, savename)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    
    window_sizes = [8,10,46,55];
    load([kdata_folder 'SNR_measurements_slices' num2str(slices)], 'S_SILVER', 'S_GR', 'S_U', 'N_SILVER', 'N_GR', 'N_U', 'SNR_SILVER', 'SNR_GR', 'SNR_U')

    for n = 1:4
        ws = window_sizes(n);
        switch ws
            case 8
                sensmask = sensmask_lowres;
                labl = '(Fibonacci)';
            case 10
                sensmask = sensmask_lowres;
                labl = [];
            case 46
                sensmask = sensmask_highres;
                labl = [];
            case 55
                sensmask = sensmask_highres;
                labl = '(Fibonacci)';
        end

        N_GR_sim = [];
        N_S_sim = [];
        N_U_sim = [];   

        for m = 1:length(slices)
            sl = slices(m);
            load([kdata_folder 'Noise_Recons_sl' num2str(sl) '_ws' num2str(ws) '.mat'], 'recon_n_GR', 'recon_n_S','recon_n_U')
            N_GR_sim = cat(2,N_GR_sim,std(recon_n_GR,[],4));
            N_S_sim = cat(2,N_S_sim,std(recon_n_S,[],4));
            N_U_sim = cat(2,N_U_sim,std(recon_n_U,[],4));
        end
  
        subplot(length(window_sizes),1,n)
        imagesc(cat(1,cat(2,sensmask.*N_U{n}./mean(N_U{n}(sensmask)), sensmask.*N_GR{n}./mean(N_U{n}(sensmask)), sensmask.*N_SILVER{n}./mean(N_U{n}(sensmask))),...
                cat(2,sensmask.*N_U_sim./mean(N_U_sim(sensmask)), sensmask.*N_GR_sim./mean(N_U_sim(sensmask)), sensmask.*N_S_sim./mean(N_U_sim(sensmask)))))
        axis image
        axis off
        rectangle('Position',[1,1,size(sensmask,2),2*size(sensmask,1)],'edgecolor', 'w','linewidth',3)
        rectangle('Position',[1+size(sensmask,2),1,size(sensmask,2),2*size(sensmask,1)],'edgecolor', 'w','linewidth',3)
        rectangle('Position',[1+size(sensmask,2)*2,1,size(sensmask,2),2*size(sensmask,1)],'edgecolor', 'w', 'linewidth',3')
        title(['N = ' num2str(window_sizes(n)) ' ' labl], 'FontSize', 12)
        text(-5, size(sensmask,1)/2, 'Measured', 'HorizontalAlignment', 'right','Color', [0  0 1 0.5])
        text(-5, size(sensmask,1)*1.5, 'Predicted', 'HorizontalAlignment', 'right','Color', [0 0 0 0.5])
        
        text(size(sensmask,2)*0.5,size(sensmask,1), 'Uniform', 'HorizontalAlignment', 'center','color','w')
        text(size(sensmask,2)*1.5,size(sensmask,1), 'GR', 'HorizontalAlignment', 'center','color','w')
        text(size(sensmask,2)*2.5,size(sensmask,1), 'SILVER', 'HorizontalAlignment', 'center','color','w')
        
        caxis([0,3])
        set(gcf, 'Position', [514 1 553 804])
        clear('GR_n_pred', 'SILVER_n_pred', 'U_n_pred')
        colormap hot
        
    end
    cb = colorbar;
    cb.Position = [0.75 0.33 0.03 0.33];

    saveas(gcf, [ savename '.svg'])
    savefig(gcf, [savename '.fig'])
end


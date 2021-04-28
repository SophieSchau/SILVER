function [] = SNR_quali_fig(sensmask_lowres, sensmask_highres, slices, kdata_folder, savename)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    
    window_sizes = [8,10,46,55];
    load([kdata_folder 'SNR_measurements_slices' num2str(slices)], 'N_SILVER', 'N_GR', 'N_U')

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
            N_GR_sim = cat(2,N_GR_sim,imrotate(std(recon_n_GR,[],4),-90));
            N_S_sim = cat(2,N_S_sim,imrotate(std(recon_n_S,[],4),-90));
            N_U_sim = cat(2,N_U_sim,imrotate(std(recon_n_U,[],4),-90));
        end
        
        sensmask = reshape(imrotate(reshape(sensmask, size(sensmask,1),size(sensmask,1),[]),-90),size(sensmask,1),[]);

        N_U{n} = reshape(imrotate(reshape(N_U{n}, size(sensmask,1),size(sensmask,1),[]),-90),size(sensmask,1),[]);
        N_GR{n} = reshape(imrotate(reshape(N_GR{n}, size(sensmask,1),size(sensmask,1),[]),-90),size(sensmask,1),[]);
        N_SILVER{n} = reshape(imrotate(reshape(N_SILVER{n}, size(sensmask,1),size(sensmask,1),[]),-90),size(sensmask,1),[]);
        
        
        
        subplot(length(window_sizes),1,n)
        imagesc(cat(1,cat(2,sensmask.*N_U{n}./mean(N_U{n}(sensmask)), sensmask.*N_GR{n}./mean(N_U{n}(sensmask)), sensmask.*N_SILVER{n}./mean(N_U{n}(sensmask))),...
                cat(2,sensmask.*N_U_sim./mean(N_U_sim(sensmask)), sensmask.*N_GR_sim./mean(N_U_sim(sensmask)), sensmask.*N_S_sim./mean(N_U_sim(sensmask)))))
        axis image
        axis off
        rectangle('Position',[1,1,size(sensmask,2),2*size(sensmask,1)],'edgecolor', 'w','linewidth',3)
        rectangle('Position',[1+size(sensmask,2),1,size(sensmask,2),2*size(sensmask,1)],'edgecolor', 'w','linewidth',3)
        rectangle('Position',[1+size(sensmask,2)*2,1,size(sensmask,2),2*size(sensmask,1)],'edgecolor', 'w', 'linewidth',3')
        title(['N = ' num2str(window_sizes(n)) ' ' labl], 'FontSize', 12)
        text(-5, size(sensmask,1)/2, 'Measured', 'HorizontalAlignment', 'right')
        text(-5, size(sensmask,1)*1.5, 'Predicted', 'HorizontalAlignment', 'right')
        
        text(size(sensmask,2)*0.5,size(sensmask,1), 'Uniform', 'HorizontalAlignment', 'center','color','w')
        text(size(sensmask,2)*1.5,size(sensmask,1), 'GR', 'HorizontalAlignment', 'center','color','w')
        text(size(sensmask,2)*2.5,size(sensmask,1), 'SILVER', 'HorizontalAlignment', 'center','color','w')
        
        caxis([0,3])
        set(gcf, 'Position', [514 1 317 804])
        clear('GR_n_pred', 'SILVER_n_pred', 'U_n_pred')
        colormap hot
        
    end
    cb = colorbar;
    cb.Position = [0.9 0.33 0.03 0.33];
    set(gcf, 'InvertHardcopy', 'off')
    set(gcf, 'Color', 'w');
    saveas(gcf, [ savename '.svg'])
    saveas(gcf, [ savename '.tiff'])
    savefig(gcf, [savename '.fig'])
end


function [] = SNR_quanti_fig(sensmask_lowres, sensmask_highres, slices, kdata_folder, savename)
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
            N_GR_sim = cat(2,N_GR_sim,std(recon_n_GR,[],4));
            N_S_sim = cat(2,N_S_sim,std(recon_n_S,[],4));
            N_U_sim = cat(2,N_U_sim,std(recon_n_U,[],4));
        end
  
        
        % significance testing (paired)
        [U_GR_p, U_GR_h] = signrank(N_GR{n}(sensmask), N_U{n}(sensmask), 'alpha', 0.05/12);
        [U_S_p, U_S_h] = signrank(N_U{n}(sensmask), N_SILVER{n}(sensmask), 'alpha', 0.05/12);
        [GR_S_p, GR_S_h] = signrank(N_GR{n}(sensmask), N_SILVER{n}(sensmask), 'alpha', 0.05/12);

        
        subplot(length(window_sizes),4,(n*4-3):(n*4-1))
        [h1,L1] = violin([N_U{n}(sensmask)./median(N_U{n}(sensmask)) N_GR{n}(sensmask)./median(N_U{n}(sensmask)) N_SILVER{n}(sensmask)./median(N_U{n}(sensmask)) ], 'mc', [], 'medc', 'k', 'facecolor', [0.1 0.1 0.1;0.1 0.1 0.1;0.1 0.1 0.1], 'side', 'left', 'bw', 0.1);
        hold on
        [h2,L2] = violin([N_U_sim(sensmask)./median(N_U_sim(sensmask)) N_GR_sim(sensmask)./median(N_U_sim(sensmask)) N_S_sim(sensmask)./median(N_U_sim(sensmask)) ], 'mc', [], 'medc', 'k', 'facecolor', [0.9 0.9 0.9;0.9 0.9 0.9;0.9 0.9 0.9], 'side', 'right','bw', 0.1);
        
        if ~U_GR_h
            line([1,2], [3.3 3.3], 'color', 'k', 'handlevisibility', 'off')
            text(1.5, 3.3, 'n.s.', 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        end
        if ~U_S_h
            line([1,3], [3.4 3.4], 'color', 'k', 'handlevisibility', 'off')
            text(2, 3.4, 'n.s.', 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        end
        if ~GR_S_h
            line([2,3], [3.35 3.35], 'color', 'k', 'handlevisibility', 'off')
            text(2.5, 3.35, 'n.s.', 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        end
        
        title(['N = ' num2str(window_sizes(n))])
        xticks(1:3)
        xticklabels({'Uniform', 'GR', 'SILVER'})
%          [h1,L1] = violin([N_U{n}(sensmask)./median(N_U{n}(sensmask)), N_U_sim(sensmask)./median(N_U_sim(sensmask))], 'mc', [], 'medc', 'b', 'facecolor', [0 0.5 1], 'facealpha', 0.3);
%          [h2,L2] = violin([N_GR{n}(sensmask)./median(N_U{n}(sensmask)),N_GR_sim(sensmask)./median(N_U_sim(sensmask))], 'mc', [], 'medc', 'r', 'facecolor', [1 0.5 0], 'facealpha', 0.3);
%          [h3,L3] = violin([N_SILVER{n}(sensmask)./median(N_U{n}(sensmask)),N_S_sim(sensmask)./median(N_U_sim(sensmask))], 'mc', [], 'medc', 'k', 'facecolor', [0.5 0.5 0.5], 'facealpha', 0.3);
         h3 = plot(xlim, [median(N_U{n}(sensmask)./median(N_U{n}(sensmask))) median(N_U{n}(sensmask)./median(N_U{n}(sensmask)))], 'linestyle', ':', 'color', [0 0.5 1], 'linewidth', 2,'DisplayName', 'Uniform (measured) median');
         h4 = plot(xlim, [median(N_GR{n}(sensmask)./median(N_U{n}(sensmask))) median(N_GR{n}(sensmask)./median(N_U{n}(sensmask)))], 'linestyle', ':', 'color', [1 0.5 0], 'linewidth', 2,'DisplayName', 'GR (measured) median');
         h5 = plot(xlim, [median(N_SILVER{n}(sensmask)./median(N_U{n}(sensmask))) median(N_SILVER{n}(sensmask)./median(N_U{n}(sensmask)))], 'linestyle', ':', 'color', [0.5 0.5 0.5], 'linewidth', 2,'DisplayName', 'SILVER (measured) median');
            
         text(3.5, 2, ['GR: ' num2str((median(N_GR{n}(sensmask)./median(N_U{n}(sensmask))))*100,3) '%'], 'Color', [1 0.5 0])
         text(3.5, 1.5, ['SILVER: ' num2str((median(N_SILVER{n}(sensmask)./median(N_U{n}(sensmask))))*100,3) '%'], 'color', [0.5 0.5 0.5])
         text(3.5, 1, ['UNIFORM: ' num2str((median(N_U{n}(sensmask)./median(N_U{n}(sensmask))))*100,3) '%'], 'color', [0 0.5 1])
         
%          legend([h1(1), h2(1),  h3(1)],{'UNIFORM', 'GR', 'SILVER'}, 'Location', 'eastoutside');
%         xticks(1:2)
%         xticklabels({'measured noise', 'predicted noise'})
        title(['N = ' num2str(window_sizes(n)) ' ' labl], 'FontSize', 12)
        ylabel('noise, normalized', 'FontSize',12)
        axis([xlim*1.15 0 3.5]) 
        set(gcf, 'Position', [1 1 904 804])
        legend('off')
    end
    subplot(length(window_sizes),4,4)
    ax = gca;
    axis off
    
    ll = legend([h1(1),   h2(1), h3, h4, h5], 'Location', 'eastoutside', 'EdgeColor', 'none');
    ll.String{1} = 'Measured noise';
    ll.String{2} = 'Predicted noise';
    ll.Position = get(ax,'position');
    saveas(gcf, [ savename '.svg'])
    saveas(gcf, [ savename '.tiff'])
    savefig(gcf, [savename '.fig'])
end


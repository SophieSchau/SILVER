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
  
        subplot(length(window_sizes)+1,1,n)
        [h1,L1] = violin([N_U{n}(sensmask)./median(N_U{n}(sensmask)) N_GR{n}(sensmask)./median(N_U{n}(sensmask)) N_SILVER{n}(sensmask)./median(N_U{n}(sensmask)) ], 'mc', [], 'medc', 'b', 'facecolor', 'b', 'side', 'left');
        hold on
        [h2,L2] = violin([N_U_sim(sensmask)./median(N_U_sim(sensmask)) N_GR_sim(sensmask)./median(N_U_sim(sensmask)) N_S_sim(sensmask)./median(N_U_sim(sensmask)) ], 'mc', [], 'medc', 'k', 'facecolor', 'k', 'side', 'right');
        title(['N = ' num2str(window_sizes(n))])
        xticks(1:3)
        xticklabels({'Uniform', 'GR', 'SILVER'})
%          [h1,L1] = violin([N_U{n}(sensmask)./median(N_U{n}(sensmask)), N_U_sim(sensmask)./median(N_U_sim(sensmask))], 'mc', [], 'medc', 'b', 'facecolor', [0 0.5 1], 'facealpha', 0.3);
%          [h2,L2] = violin([N_GR{n}(sensmask)./median(N_U{n}(sensmask)),N_GR_sim(sensmask)./median(N_U_sim(sensmask))], 'mc', [], 'medc', 'r', 'facecolor', [1 0.5 0], 'facealpha', 0.3);
%          [h3,L3] = violin([N_SILVER{n}(sensmask)./median(N_U{n}(sensmask)),N_S_sim(sensmask)./median(N_U_sim(sensmask))], 'mc', [], 'medc', 'k', 'facecolor', [0.5 0.5 0.5], 'facealpha', 0.3);
         h3 = plot(xlim, [median(N_U{n}(sensmask)./median(N_U{n}(sensmask))) median(N_U{n}(sensmask)./median(N_U{n}(sensmask)))], 'linestyle', ':', 'color', [0 0.5 1], 'linewidth', 1,'DisplayName', 'Uniform (measured) median');
         h4 = plot(xlim, [median(N_GR{n}(sensmask)./median(N_U{n}(sensmask))) median(N_GR{n}(sensmask)./median(N_U{n}(sensmask)))], 'linestyle', ':', 'color', [1 0.5 0], 'linewidth', 1,'DisplayName', 'GR (measured) median');
         h5 = plot(xlim, [median(N_SILVER{n}(sensmask)./median(N_U{n}(sensmask))) median(N_SILVER{n}(sensmask)./median(N_U{n}(sensmask)))], 'linestyle', ':', 'color', [0.5 0.5 0.5], 'linewidth', 1,'DisplayName', 'SILVER (measured) median');
%          legend([h1(1), h2(1),  h3(1)],{'UNIFORM', 'GR', 'SILVER'}, 'Location', 'eastoutside');
%         xticks(1:2)
%         xticklabels({'measured noise', 'predicted noise'})
        title(['N = ' num2str(window_sizes(n)) ' ' labl], 'FontSize', 12)
        ylabel('noise, normalized', 'FontSize',12)
        axis([xlim 0 3]) 
        set(gcf, 'Position', [1 1 512 804])
        legend('off')
    end
    subplot(length(window_sizes)+1,1,n+1)
    ax = gca;
    axis off
    
    ll = legend([h1(1),   h2(1), h3, h4, h5], 'Location', 'eastoutside', 'EdgeColor', 'none');
    ll.String{1} = 'Measured noise';
    ll.String{2} = 'Predicted noise';
    ll.Position = get(ax,'position');
    saveas(gcf, [ savename '.svg'])
    savefig(gcf, [savename '.fig'])
end


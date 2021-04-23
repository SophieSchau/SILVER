function [efficiency,noise, alpha] = efficiency_and_noise_amplification(NSpokes,efficiency_metric, mat_size, savefilename, alpha)
%EFFICIENCY_AND_NOISE_AMPLIFICATION How does uniformity relate to noise amplification?
%   For a fixed number of spokes (NSpokes), how does the efficiency
%   (calculated with some efficiency metric) correlate with noise
%   amplification when reconstructed into an image of matrix size mat_size
%   x mat_size using pseudo-inversion  + a final masking step to zero out 
%   the un-measured corners of k-space. 100 iterations of noise data are
%   calculated and the pixel wise standard deviation is used as a noise
%   estimate.
%
%   The variable savefilename can be used to save the results to a file
%   with the defined path and name.

    if nargin < 4
        savefilename = [];
    end

    if nargin < 5
        rng(1)
        alpha = [gr2D, rand(1,99)];
    end

    for i = 1:mat_size
        for j = 1:mat_size
            mask(i,j) = (i-(floor(mat_size/2)+1))^2+(j-(floor(mat_size/2)+1))^2-floor(mat_size/2)^2<0;
        end
    end


    if exist([savefilename '.mat'], 'file')
        load(savefilename)
        warning('using previously calculated values')
    else




        for n = 1:length(alpha)
            efficiency(n) = efficiency_2D([0:NSpokes-1].*pi*alpha(n), efficiency_metric);

            k = reshape(gen_radial(alpha(n)*180,mat_size*2,NSpokes,1,360,1),[],1,2);
            k = repmat(k,1,100,1);
            E = xfm_NUFFT([mat_size,mat_size,1,100],ones(mat_size),[],k, 'wi', 1);
            nn = (randn(size(k,1),100)+1i*randn(size(k,1),100))./sqrt(2);
            recon_n = reshape(E.iter(nn,@pcg,0,100), mat_size, mat_size, 1, 100);

            for i = 1:100
                recon_n(:,:,:,i) = ifft2(ifftshift(fftshift(fft2(recon_n(:,:,:,i))).*mask));
            end

            noise(:,:,n) = std(recon_n,[],4);

        end
    end


    [noise_means, sortorder] = sort(squeeze(mean(noise,[1,2])));
    efficiency = efficiency(sortorder);
    alpha = alpha(sortorder);
    noise = noise(:,:,sortorder);



    idx = true(size(alpha));%efficiency>0.1;


    figure
    scatter(efficiency(idx),noise_means(idx), 30,alpha(idx),'filled')
    hold on


    [f,gof] = fit(efficiency(idx)',noise_means(idx),'poly1');
    h = plot(f);
    set(h, 'color', 'k')
    legend(h,['R^2 = ' num2str(gof.rsquare)], 'box', 'off')
    h2 = colorbar;
    ylabel(h2, '\alpha')
    xlabel('efficiency, \eta')
    ylabel('noise')
    set(gca, 'FontSize', 14)
    title([num2str(NSpokes) ' spokes - ' num2str(mat_size) 'x' num2str(mat_size) ' matrix size'])


    if ~isempty(savefilename)
        save(savefilename,'efficiency','noise', 'alpha')


    end
end

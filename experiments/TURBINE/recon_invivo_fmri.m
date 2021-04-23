function [] = recon_invivo_fmri(slices,sens_lowres, sens_highres, kdata_folder)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here


    experiments = [1:2];
    exp_names = {'GR_rs', 'SILVER_rs'};
    ratios = [...
                gr2D,...
                SILVER_2D([10,46], 'electrostatic_potential')];
    window_sizes_per_exp = {[10,46],[10,46]};
   %%    
    for ex = experiments
        ratio = ratios(ex);
        window_sizes = window_sizes_per_exp{ex};


        for sl = slices
            kdata_file_a = dir([kdata_folder 'meas_*_' exp_names{ex} '*_slice' num2str(sl) '.mat']);
            load([kdata_folder  kdata_file_a.name], 'kdata')
            kdata_a = kdata;    

            for ws = window_sizes
                Nframes = floor(size(kdata_a,2)/ws);

                kspace = gen_radial(roundn(ratio*180,-3),200, ws*Nframes,1,360,1);
                kdata_rs_a = kdata_a(:,1:ws*Nframes,:,:);


                if ws == 8 || ws == 10
                    sens = sens_lowres;
                    kspace = kspace./0.3;
                    idx = find(abs(kspace(:,1,1))<pi);
                    kspace = reshape(kspace(idx,:,:), [], Nframes,2);
                    kdata_rs_a = reshape(kdata_rs_a(idx,:,:,:),[],Nframes, 32);
                    mat_size = 30;

                else
                    sens = sens_highres;
                    idx = find(abs(kspace(:,1,1))<pi);
                    kspace = reshape(kspace(idx,:,:), [], Nframes,2);
                    kdata_rs_a = reshape(kdata_rs_a(idx,:,:,:),[],Nframes, 32);
                    mat_size = 100;
                end

                mask = zeros(mat_size);
                for i = 1:mat_size
                    for j = 1:mat_size
                        mask(i,j) = (i-(floor(mat_size/2)+1))^2+(j-(floor(mat_size/2)+1))^2-floor(mat_size/2)^2<0;
                    end
                end



                E = xfm_NUFFT([mat_size,mat_size,1,Nframes],sens(:,:,sl,:),[],kspace, 'wi', 1);
                recon_a = reshape(E.iter(kdata_rs_a,@pcg,0,100), mat_size, mat_size, 1, []);   

                for t = 1:size(recon_a,4)
                    recon_a(:,:,:,t) = ifft2(ifftshift(fftshift(fft2(recon_a(:,:,:,t))).*mask));
                end

                save([kdata_folder 'recons/Recon_sl' num2str(sl) '_ws' num2str(ws) kdata_file_a.name(1:end-12) '_iter.mat'], 'recon_a')
            end
        end
    end






end


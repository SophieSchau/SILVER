function [] = recon_invivo(slices,sens_lowres, sens_highres, kdata_folder)
    experiments = [1:7];
    exp_names = {'UNI*8', 'UNI*10', 'UNI*46', 'UNI*55', 'GR', 'SILVER*8', 'SILVER*10'};
    ratios = [...
                1/8,...
                1/10,...
                1/46,...
                1/55,...
                gr2D,...
                SILVER_2D([8,55], 'electrostatic_potential'),...
                SILVER_2D([10,46], 'electrostatic_potential')];

    window_sizes_per_exp = {8, 10, 46, 55, [8,10,46,55],[8,55],[10,46]};






    %%    
    for ex = experiments
        ratio = ratios(ex);
        window_sizes = window_sizes_per_exp{ex};


        for sl = slices
            kdata_file_a = dir([kdata_folder 'meas_*_' exp_names{ex} '*a_slice' num2str(sl) '.mat']);
            load([kdata_folder  kdata_file_a.name], 'kdata')
            kdata_a = kdata;

            kdata_file_b = dir([kdata_folder 'meas_*_' exp_names{ex} '*b_slice' num2str(sl) '.mat']);
            load([kdata_folder  kdata_file_b.name], 'kdata')
            kdata_b = kdata;
            
            for spoke = 1:size(kdata_a,2)
                ph = angle(mean(conj(kdata_a(:,spoke,:,:)).*kdata_b(:,spoke,:,:), 'all'));
                kdata_a(:,spoke,:,:) = kdata_a(:,spoke,:,:)*exp(1j*ph);
            end     

            for ws = window_sizes
                if exist([kdata_folder 'recons/Recon_sl' num2str(sl) '_ws' num2str(ws) kdata_file_a.name(1:end-12) '_iter.mat'], 'file')
                    warning(['A file named "' kdata_folder 'recons/Recon_sl' num2str(sl) '_ws' num2str(ws) kdata_file_a.name(1:end-12) '_iter.mat" already exists! no reconstruction performed'])
                    continue
                end
                Nframes = floor(size(kdata_a,2)/ws);

                kspace = gen_radial(roundn(ratio*180,-3),200, ws*Nframes,1,360,1);
                kdata_rs_a = kdata_a(:,1:ws*Nframes,:,:);
                kdata_rs_b = kdata_b(:,1:ws*Nframes,:,:);


                if ws == 8 || ws == 10
                    sens = sens_lowres;
                    kspace = kspace./0.3;
                    idx = find(abs(kspace(:,1,1))<pi);
                    kspace = reshape(kspace(idx,:,:), [], Nframes,2);
                    kdata_rs_a = reshape(kdata_rs_a(idx,:,:,:),[],Nframes, 32);
                    kdata_rs_b = reshape(kdata_rs_b(idx,:,:,:),[],Nframes, 32);
                    mat_size = 30;

                else
                    sens = sens_highres;
                    idx = find(abs(kspace(:,1,1))<pi);
                    kspace = reshape(kspace(idx,:,:), [], Nframes,2);
                    kdata_rs_a = reshape(kdata_rs_a(idx,:,:,:),[],Nframes, 32);
                    kdata_rs_b = reshape(kdata_rs_b(idx,:,:,:),[],Nframes, 32);
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
                recon_b = reshape(E.iter(kdata_rs_b,@pcg,0,100), mat_size, mat_size, 1, []);   

                for t = 1:size(recon_a,4)
                    recon_a(:,:,:,t) = ifft2(ifftshift(fftshift(fft2(recon_a(:,:,:,t))).*mask));
                    recon_b(:,:,:,t) = ifft2(ifftshift(fftshift(fft2(recon_b(:,:,:,t))).*mask));
                end

                save([kdata_folder 'recons/Recon_sl' num2str(sl) '_ws' num2str(ws) kdata_file_a.name(1:end-12) '_iter.mat'], 'recon_a', 'recon_b')
            end
        end
    end






end
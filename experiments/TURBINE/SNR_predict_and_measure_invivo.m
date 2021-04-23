function [] = SNR_predict_and_measure_invivo(slices, Nspokes, sens_lowres, sens_highres, ncov1, kdata_folder)
exp_names = {'UNI*8', 'UNI*10', 'UNI*46', 'UNI*55', 'GR', 'SILVER*8', 'SILVER*10'};
ratios = [...
            1/8,...
            1/10,...
            1/46,...
            1/55,...
            gr2D,...
            SILVER_2D([8,55], 'electrostatic_potential'),...
            SILVER_2D([10,46], 'electrostatic_potential')];
            

window_sizes= [8, 10, 46, 55];
 
%% measurements 
if ~exist([kdata_folder 'SNR_measurements_slices' num2str(slices) '.mat'], 'file')
    for n = 1:length(window_sizes)
        ws = window_sizes(n);
        switch ws
            case 8
                U_ex = 'UNI_8a';
                SILVER_ex = 'SILVER_8*a';
            case 10
                U_ex = 'UNI_10a';
                SILVER_ex = 'SILVER_10*a';
            case 46
                U_ex = 'UNI_46a';
                SILVER_ex = 'SILVER_10*a';
            case 55
                U_ex = 'UNI_55a';
                SILVER_ex = 'SILVER_8*a';
        end
        
        
        GR_a = [];
        GR_b = [];
        U_a = [];
        U_b = [];
        SILVER_a = [];
        SILVER_b = [];
        
        for m = 1:length(slices)
            kdata_file_GR = dir([kdata_folder '../meas_*GRa*slice' num2str(slices(m)) '.mat']);
            load([kdata_folder 'Recon_sl' num2str(slices(m)) '_ws' num2str(ws) kdata_file_GR.name(1:end-12) '_iter.mat'], 'recon_a', 'recon_b')
            GR_a = cat(2,GR_a,recon_a(:,:,:,1:end));
            GR_b = cat(2,GR_b,recon_b(:,:,:,1:end));

            kdata_file_U = dir([kdata_folder '../meas_*_' U_ex '*slice' num2str(slices(m)) '.mat']);
            load([kdata_folder 'Recon_sl' num2str(slices(m)) '_ws' num2str(ws) kdata_file_U.name(1:end-12) '_iter.mat'], 'recon_a', 'recon_b')
            U_a = cat(2,U_a,recon_a(:,:,:,1:end));
            U_b = cat(2,U_b,recon_b(:,:,:,1:end));

            kdata_file_SILVER = dir([kdata_folder '../meas_*_' SILVER_ex '*slice' num2str(slices(m)) '.mat']);
            load([kdata_folder 'Recon_sl' num2str(slices(m)) '_ws' num2str(ws) kdata_file_SILVER.name(1:end-12) '_iter.mat'], 'recon_a', 'recon_b')
            SILVER_a = cat(2,SILVER_a,recon_a(:,:,:,1:end));
            SILVER_b = cat(2,SILVER_b,recon_b(:,:,:,1:end));

        end
        
        figure
        imagesc(cat(1,abs(mean(cat(1,cat(2,GR_a,GR_b),cat(2,SILVER_a,SILVER_b),cat(2,U_a,U_b)),4)),...
            abs(mean(cat(1,cat(2,(GR_a+GR_b)./2,GR_a-GR_b),cat(2,(SILVER_a+SILVER_b)./2,SILVER_a-SILVER_b),cat(2,(U_a+U_b)./2,U_a-U_b)),4))))
        drawnow
        S_GR{n} = mean(abs(GR_a)+abs(GR_b),4);
        N_GR{n} = std(abs(GR_a)-abs(GR_b),[],4);
        SNR_GR{n} = S_GR{n}./(sqrt(2)*N_GR{n});

        S_U{n} = mean(abs(U_a)+abs(U_b),4);
        N_U{n} = std(abs(U_a)-abs(U_b),[],4);
        SNR_U{n} = S_U{n}./(sqrt(2)*N_U{n});

        S_SILVER{n} = mean(abs(SILVER_a)+abs(SILVER_b),4);
        N_SILVER{n} = std(abs(SILVER_a)-abs(SILVER_b),[],4);
        SNR_SILVER{n} = S_SILVER{n}./(sqrt(2)*N_SILVER{n});
        
    end
    save([kdata_folder 'SNR_measurements_slices' num2str(slices)], 'S_SILVER', 'S_GR', 'S_U', 'N_SILVER', 'N_GR', 'N_U', 'SNR_SILVER', 'SNR_GR', 'SNR_U')

else
    load([kdata_folder 'SNR_measurements_slices' num2str(slices)], 'S_SILVER', 'S_GR', 'S_U', 'N_SILVER', 'N_GR', 'N_U', 'SNR_SILVER', 'SNR_GR', 'SNR_U')
    warning('using pre-measured SNR values')
end


%% prediction - monte carlo

for n = 1:length(window_sizes)
    ws = window_sizes(n);

    Nframes = floor(Nspokes/ws);
    
    if ws == 8 || ws == 55
        SILVER_ex = 6;
    else
        SILVER_ex = 7;
    end
    
    kspace_GR = gen_radial(roundn(gr2D*180,-3),200, ws*Nframes,1,360,1);
    kspace_S = gen_radial(roundn(ratios(SILVER_ex)*180,-3),200, ws*Nframes,1,360,1);
    kspace_U = gen_radial(roundn(1/ws*180,-3),200, ws*Nframes,1,360,1);



    if ws == 8 || ws == 10
        sens = sens_lowres;
        mat_size = 30;
        sensmask{n} = logical(abs(mean(reshape(sens(:,:,slices,:),mat_size, mat_size*length(slices),1,[]).^2,4)));

        kspace_GR = kspace_GR./0.3;
        kspace_S = kspace_S./0.3;
        kspace_U = kspace_U./0.3;
        idx = find(abs(kspace_GR(:,1,1))<pi);
        kspace_GR = reshape(kspace_GR(idx,:,:), [], Nframes,2);
        kspace_S = reshape(kspace_S(idx,:,:), [], Nframes,2);
        kspace_U = reshape(kspace_U(idx,:,:), [], Nframes,2);

        mat_size = 30;

    else
        sens = sens_highres;
        mat_size = 100;
        sensmask{n} = logical(abs(mean(reshape(sens(:,:,slices,:),mat_size, mat_size*length(slices),1,[]).^2,4)));

        idx = find(abs(kspace_GR(:,1,1))<pi);
        kspace_GR = reshape(kspace_GR(idx,:,:), [], Nframes,2);
        kspace_S = reshape(kspace_S(idx,:,:), [], Nframes,2);
        kspace_U = reshape(kspace_U(idx,:,:), [], Nframes,2);
        mat_size = 100;
    end
            
    mask = zeros(mat_size);
    for i = 1:mat_size
        for j = 1:mat_size
            mask(i,j) = (i-(floor(mat_size/2)+1))^2+(j-(floor(mat_size/2)+1))^2-floor(mat_size/2)^2<0;
        end
    end



    for sl  = slices
        if ~exist([kdata_folder 'Noise_Recons_sl' num2str(sl) '_ws' num2str(ws) '.mat'], 'file')

            E_GR = xfm_NUFFT([mat_size,mat_size,1,Nframes],sens(:,:,sl,:),[],kspace_GR, 'wi', 1);
            E_S = xfm_NUFFT([mat_size,mat_size,1,Nframes],sens(:,:,sl,:),[],kspace_S, 'wi', 1);
            E_U = xfm_NUFFT([mat_size,mat_size,1,Nframes],sens(:,:,sl,:),[],kspace_U, 'wi', 1);


            nn = (randn(length(idx)*ws*Nframes,32)+1i*randn(length(idx)*ws*Nframes,32))./sqrt(2);
            c = sqrtm(ncov1);
            nc = nn*c;
            kdata_n = reshape(nc,[],Nframes,32);


            recon_n_GR = reshape(E_GR.iter(kdata_n,@pcg,0,100), mat_size, mat_size, 1, Nframes);
            recon_n_S = reshape(E_S.iter(kdata_n,@pcg,0,100), mat_size, mat_size, 1, Nframes);
            recon_n_U = reshape(E_U.iter(kdata_n,@pcg,0,100), mat_size, mat_size, 1, Nframes);



            for n = 1:Nframes
                recon_n_GR(:,:,:,n) = ifft2(ifftshift(fftshift(fft2(recon_n_GR(:,:,:,n))).*mask));
                recon_n_S(:,:,:,n) = ifft2(ifftshift(fftshift(fft2(recon_n_S(:,:,:,n))).*mask));
                recon_n_U(:,:,:,n) = ifft2(ifftshift(fftshift(fft2(recon_n_U(:,:,:,n))).*mask));

            end

            save([kdata_folder 'Noise_Recons_sl' num2str(sl) '_ws' num2str(ws) '.mat'], 'recon_n_GR', 'recon_n_S','recon_n_U')
        end
    end
end

end
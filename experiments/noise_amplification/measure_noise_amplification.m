function [] = measure_noise_amplification(sens, ncov, N, S, efficiency_metric, start_angle, savename)
%MEASURE_NOISE_AMPLIFICATION How does the trajectory (SILVER, golden ratio, or 
%   uniform radial) affect the noise amplification maps for a set 
%   sensitivity profile?
%                                              
%   Sophie Schauman 2020                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_size = [size(sens,1) size(sens,2) size(sens,3)];
Nframes = 100;
% Do SILVER optimization for the sets of ranges
for s = S
   savename_ratio = ['experiments/precalculated_' efficiency_metric '/silver_' strrep(num2str(s{:}),' ', '_') '.mat'];
   if ~exist(savename_ratio, 'file')
        SILVER_2D(s{:},efficiency_metric, savename_ratio) ;
   end
   disp([num2str(s{:}) ' done'])
end

% Generate g-factor and noise amplification maps

for i = 1:mat_size(1)
    for j = 1:mat_size(2)
        mask(i,j) = (i-(floor(mat_size(1)/2)+1))^2+(j-(floor(mat_size(1)/2)+1))^2-floor(mat_size(1)/2)^2<0;
    end
end

rng(start_angle+1)
nn = (randn(200*N*Nframes,size(sens,4))+1i*randn(200*N*Nframes,size(sens,4)))./sqrt(2);
c = sqrtm(ncov);
nc = nn*c;
nnn = reshape(nc,[],Nframes,size(sens,4));
            
for m = 1:length(S)
    for n = 1:length(S{m})
        if S{m}(n) ==  N

            if ~exist([savename 'SILVER_N_' num2str(N) '_S_' strrep(num2str(S{m}),' ', '_') '.mat'], 'file')
                load(['experiments/precalculated_' efficiency_metric '/silver_' strrep(num2str(S{m}),' ', '_') '.mat'],'ratio');
                k = reshape(gen_radial(ratio*180,200,N,1,360,1),[],2);
                rotmat = [cosd(start_angle) -sind(start_angle);sind(start_angle) cosd(start_angle)];
                k = k*rotmat;
                k = reshape(k,[],1,2);
                k = repmat(k,1,Nframes,1);

                E = xfm_NUFFT([mat_size,Nframes],sens,[],k, 'wi', 1);
                recon_n = reshape(E.iter(nnn,@pcg,0,Nframes), mat_size(1), mat_size(2), mat_size(3), Nframes);

                for i = 1:Nframes
                    recon_n(:,:,:,i) = ifft2(ifftshift(fftshift(fft2(recon_n(:,:,:,i))).*mask));
                end

                SILVER_n = std(recon_n,[],4);
                
                save([savename 'SILVER_N_' num2str(N) '_S_' strrep(num2str(S{m}),' ', '_') '.mat'],'SILVER_n','recon_n' )
            end
        end
    end
end


if ~exist([savename 'GR_N_' num2str(N) '.mat'], 'file')
    ratio = gr2D;
    k = reshape(gen_radial(ratio*180,200,N,1,360,1),[],2);
    rotmat = [cosd(start_angle) -sind(start_angle);sind(start_angle) cosd(start_angle)];
    k = k*rotmat;
    k = reshape(k,[],1,2);
    k = repmat(k,1,Nframes,1);

    E = xfm_NUFFT([mat_size,Nframes],sens,[],k, 'wi', 1);
    recon_n = reshape(E.iter(nnn,@pcg,0,Nframes), mat_size(1), mat_size(2), mat_size(3), Nframes);

    for i = 1:Nframes
        recon_n(:,:,:,i) = ifft2(ifftshift(fftshift(fft2(recon_n(:,:,:,i))).*mask));
    end

    GR_n = std(recon_n,[],4);
    save([savename 'GR_N_' num2str(N) '.mat'],'GR_n','recon_n')
end

if ~exist([savename 'U_N_' num2str(N) '.mat'], 'file')
    ratio = 1/N;
    k = reshape(gen_radial(ratio*180,200,N,1,360,1),[],2);
    rotmat = [cosd(start_angle) -sind(start_angle);sind(start_angle) cosd(start_angle)];
    k = k*rotmat;
    k = reshape(k,[],1,2);
    k = repmat(k,1,Nframes,1);

    E = xfm_NUFFT([mat_size,Nframes],sens,[],k, 'wi', 1);
    recon_n = reshape(E.iter(nnn,@pcg,0,Nframes), mat_size(1), mat_size(2), mat_size(3), Nframes);

    for i = 1:Nframes
        recon_n(:,:,:,i) = ifft2(ifftshift(fftshift(fft2(recon_n(:,:,:,i))).*mask));
    end

    U_n = std(recon_n,[],4);
    save([savename 'U_N_' num2str(N) '.mat'],'U_n','recon_n')
end
       



end
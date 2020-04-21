%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   How does the trajectory (SILVER,golden ratio, or uniform radial) affect
%   the g-factor maps for a set sensitivity profile?
%                                              
%   Sophie Schauman 2020                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% 1. Choose  set of window sizes, S, to consider
S = [16,32,48];

%% 2. Do SILVER optimization for that range

savename = ['examples/precalculated/silver_' strrep(num2str(S),' ', '_') '.mat'];
if ~exist(savename, 'file')
    ratio = SILVER_2D(S,'electrostatic_potential',savename);
else
    load(savename,'ratio')
end

%% 3. Generate g-factors for SILVER, GR, and uniform sampling
savename = 'examples/example5_gfactor/example5_gfactor_result.mat';

if ~exist(savename, 'file')
    load('examples/example4_psf/compressed_sensitivities.mat', 'sens')
    psens = randn(64,64,1,8)+1i*randn(64,64,1,8);
    psens = psens./sum(abs(psens),4);
    for i = 1:8
        psens(:,:,1,i) = imresize(imrotate(sens(:,:,1,i),0),[64,64]);
    end

    for n = 1:length(S)
        k_Uniform= reshape(sqrt(2)*gen_radial_traj((0:S(n)-1)*1/S(n)*pi, 128, [])./pi,[],2);
        
        k_GR= reshape(sqrt(2)*gen_radial_traj((0:S(n)-1)*gr2D*pi, 128, [])./pi,[],2);
        
        k_SILVER= reshape(sqrt(2)*gen_radial_traj((0:S(n)-1)*ratio*pi, 128, [])./pi,[],2);
        
%         E_Uniform{n} = xfm_NUFFT([64 64 1 1],psens,[],k_Uniform);
%         E_GR{n} = xfm_NUFFT([64 64 1 1],psens,[],k_GR);
%         E_SILVER{n} = xfm_NUFFT([64 64 1 1],psens,[],k_SILVER);
%         for m = 1:1000
%             rng(m)
%             kdata = (randn(size(k_Uniform,1),size(k_Uniform,2),8) + 1i* randn(size(k_Uniform,1),size(k_Uniform,2),8))/sqrt(2);
%             recon_Uniform(:,:,m) = fista(E_Uniform{n}, kdata, 1, 0.000000, ...
%                 [64, 64, 1,1], 50, 0.5);
%             recon_GR(:,:,m) = fista(E_GR{n}, kdata, 1, 0.000000, ...
%                 [64, 64, 1,1], 50, 0.5);
%             recon_SILVER(:,:,m) = fista(E_SILVER{n}, kdata, 1, 0.000000, ...
%                 [64, 64, 1,1], 50, 0.5);
%             
%         end

        SILVER_g{n} = gFactor(psens,k_SILVER);
        GR_g{n} = gFactor(psens,k_GR);
        Uniform_g{n} = gFactor(psens,k_Uniform);
%         SILVER_g{n} = std(recon_SILVER,[],3);
%         GR_g{n} = std(recon_GR,[],3);
%         Uniform_g{n} = std(recon_Uniform,[],3);

        clear('k_Uniform', 'k_GR','k_SILVER')
    end
    try
        save(savename, 'SILVER_g', 'GR_g', 'Uniform_g')
    catch
        warning('g-factor maps not saved because the save location could not be found')
    end
else
    load(savename, 'SILVER_g', 'GR_g', 'Uniform_g')
    warning('Using precalculated g-factor maps)')
end

%% 4. Vizualise
figure(1)
[columnsInImage, rowsInImage] = meshgrid(1:64, 1:64);
centerX = 32;
centerY = 32;
radius = 60;
circleMask = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;

c = [1,100; 1,10; 1, 1.1];

for n = 1:length(S)
    subplot(length(S),3,n*3-2)
    im = abs(Uniform_g{n}).*circleMask;
    disp(['Maximum g-factor for Uniform sampling (' num2str(S(n)) ' spokes) is ' num2str(max(im(:)))])
    imagesc(im)
    axis image
    axis off
    if n ==1
        text(32,-4,'UNIFORM', 'fontsize', 20,'HorizontalAlignment','center')
    end
    caxis(c(n,:))
    colorbar
    set(gca,'FontSize', 18)
    text(-2,32, [num2str(S(n)) ' spokes'], 'fontsize', 20,'HorizontalAlignment','right')
    
    
    subplot(length(S),3,n*3-1)
    im = abs(GR_g{n}).*circleMask;
    disp(['Maximum g-factor for Golden ratio sampling (' num2str(S(n)) ' spokes) is ' num2str(max(im(:)))])
    imagesc(im)
    axis image
    axis off
    if n ==1
        text(32,-4,'GOLDEN RATIO', 'fontsize', 20,'HorizontalAlignment','center')
    end
    caxis(c(n,:))
    colorbar
    set(gca,'FontSize', 18)
    
    subplot(length(S),3,n*3)
    im = abs(SILVER_g{n}).*circleMask;
    disp(['Maximum g-factor for SILVER sampling (' num2str(S(n)) ' spokes) is ' num2str(max(im(:)))])
    imagesc(im)
    axis image
    axis off
    if n ==1
        text(32,-4,'SILVER', 'fontsize', 20,'HorizontalAlignment','center')
    end
    caxis(c(n,:));
    colorbar
    set(gca,'FontSize', 18)
    disp(' ')
    
end
% set(gcf,'Position', [440 1 794 797])
% savefig('examples/example5_gfactor/example5_gfactor_result.fig')
% saveas(gcf,'examples/example5_gfactor/example5_gfactor_result.tiff')

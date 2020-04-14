%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   How does the trajectory (SILVER,golden ratio, or uniform radial) affect
%   the g-factor maps for a set sensitivity profile?
%                                              
%   Sophie Schauman 2020                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% 1. Choose  set of window sizes, S, to consider
S = [16,32,64];

%% 2. Do SILVER optimization for that range

savename = ['examples/precalculated/silver_' strrep(num2str(S),' ', '_') '.mat'];
if ~exist(savename, 'file')
    ratio = SILVER_2D(S,'electrostatic_potential',savename);
else
    load(savename,'ratio')
end

%% 3. Generate g-factors for SILVER, GR, and uniform sampling
savename = 'examples/example5_gfactor/gfactormaps.mat';

if ~exist(savename, 'file')
    load('examples/example4_psf/compressed_sensitivities.mat','sens')
    psens = zeros(64,64,1,8);
    for i = 1:8
        psens(:,:,1,i) = imresize(sens(:,:,1,i),[64,64]);
    end

    for n = 1:length(S)
        k_Uniform= sqrt(2)*gen_radial_traj((0:S(n)-1)*1/S(n)*pi, 128, [])/(pi);
        k_GR= sqrt(2)*gen_radial_traj((0:S(n)-1)*gr2D*pi, 128, [])/(pi);
        k_SILVER= sqrt(2)*gen_radial_traj((0:S(n)-1)*ratio*pi, 128, [])/(pi);

        SILVER_g{n} = gFactor(psens,k_SILVER);
        GR_g{n} = gFactor(psens,k_GR);
        Uniform_g{n} = gFactor(psens,k_Uniform);
    end
    save(savename, 'SILVER_g', 'GR_g', 'Uniform_g')
else
    load(savename, 'SILVER_g', 'GR_g', 'Uniform_g')
    warning('Using precalculated g-factor maps)')
end

%% 4. Vizualise
figure(1)
[columnsInImage, rowsInImage] = meshgrid(1:64, 1:64);
centerX = 32;
centerY = 32;
radius = 32;
circleMask = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;

for n = 1:length(S)
    subplot(length(S),3,n*3-2)
    im = abs(Uniform_g{n}).*circleMask;
    disp(['Maximum g-factor for Uniform sampling (' num2str(S(n)) ' spokes) is ' num2str(max(im(:)))])
    imagesc(im)
    axis image
    axis off
    if n ==1
        text(32,-4,'UNIFORM', 'fontsize', 20,'HorizontalAlignment','center')
        caxis([1,200])
    elseif n == 2
        caxis([1,5])
    elseif n == 3
        caxis([1,1.4])
    end
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
        caxis([1,200])
    elseif n == 2
        caxis([1,5])
    elseif n == 3
        caxis([1,1.4])
    end
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
        caxis([1,200])
    elseif n == 2
        caxis([1,5])
    elseif n == 3
        caxis([1,1.4])
    end
    colorbar
    set(gca,'FontSize', 18)
    disp(' ')
    
end
set(gcf,'Position', [440 1 794 797])
savefig('examples/example5_gfactor/example5_gfactor_result.fig')
saveas(gcf,'examples/example5_gfactor/example5_gfactor_result.tiff')
save('examples/example5_gfactor/example5_gfactor_result.mat', 'SILVER_g', 'GR_g', 'Uniform_g')

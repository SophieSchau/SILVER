%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   What do the point-spread-functions of trajectories optimized with
%   SILVER look like compared with golden ratio method or uniform radial
%   sampling?
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
    ratio = SILVER_2D(S,'electrostatic_potential',savename) ;
else
    load(savename,'ratio')
end

%% 3. Generate PSF's for SILVER, GR, and uniform sampling

for n = 1:length(S)
    SILVER_PSF{n} = psf_radial(ratio, S(n), 64);
    GR_PSF{n} = psf_radial(gr2D, S(n), 64);
    Uniform_PSF{n} = psf_radial(1/S(n), S(n), 64);
end


%% 4. Vizualise PSFs
figure(1)
delta = zeros(64,64);
delta(32,32) = 1;

for n = 1:length(S)
    subplot(length(S),3,n*3-2)
    imagesc((abs(reshape(Uniform_PSF{n}*delta(:),64, 64))))
    axis image
    axis off
    if n ==1
        text(32,-4,'UNIFORM', 'fontsize', 20,'HorizontalAlignment','center')
        cl = caxis;
    end
    text(-2,32, [num2str(S(n)) ' spokes'], 'fontsize', 20,'HorizontalAlignment','right')
    caxis(cl)
    subplot(length(S),3,n*3-1)
    imagesc((abs(reshape(GR_PSF{n}*delta(:),64, 64))))
    axis image
    axis off
    if n ==1
        text(32,-4,'GOLDEN RATIO', 'fontsize', 20,'HorizontalAlignment','center')
    end
    caxis(cl)
    subplot(length(S),3,n*3)
    imagesc((abs(reshape(SILVER_PSF{n}*delta(:),64, 64))))
    axis image
    axis off
    if n ==1
        text(32,-4,'SILVER', 'fontsize', 20,'HorizontalAlignment','center')
    end
    caxis(cl)
end
set(gcf,'Position', [390 1 844 797])
savefig('examples/example4_psf/example4_psf_result.fig')
saveas(gcf,'examples/example4_psf/example4_psf_result.tiff')

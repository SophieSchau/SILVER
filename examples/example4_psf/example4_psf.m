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
S = [16,32,64];

%% 2. Do SILVER optimization for that range

savename = ['examples/precalculated/silver_' strrep(num2str(S),' ', '_') '.mat'];
if ~exist(savename, 'file')
    ratio = SILVER_2D(S,'electrostatic_potential',savename) ;
else
    load(savename,'ratio')
end


%% 3. Generate PSF's for SILVER, GR, and uniform sampling and vizualise

for n = 1:length(S)
    SILVER_PSF{n} = psf_radial(ratio, S(n), 64);
    GR_PSF{n} = psf_radial(gr2D, S(n), 64);
    Uniform_PSF{n} = psf_radial(1/S(n), S(n), 64);
end
    
figure(1)
for n = 1:length(S)
    subplot(length(S),3,n*3-2)
    imagesc(log(abs(Uniform_PSF{n})))
    axis image
    axis off
    if n ==1
        text(32,-4,'UNIFORM', 'fontsize', 20,'HorizontalAlignment','center')
    end
    text(-2,32, [num2str(S(n)) ' spokes'], 'fontsize', 20,'HorizontalAlignment','right')
    
    subplot(length(S),3,n*3-1)
    imagesc(log(abs(GR_PSF{n})))
    axis image
    axis off
    if n ==1
        text(32,-4,'GOLDEN RATIO', 'fontsize', 20,'HorizontalAlignment','center')
    end
    
    subplot(length(S),3,n*3)
    imagesc(log(abs(SILVER_PSF{n})))
    axis image
    axis off
    if n ==1
        text(32,-4,'SILVER', 'fontsize', 20,'HorizontalAlignment','center')
    end
end
set(gcf,'Position', [440 1 794 797])
savefig('examples/example4_psf/example4_psf_result.fig')
saveas(gcf,'examples/example4_psf/example4_psf_result.tiff')

% %% 4. Analyse PSF's
% for n = 1:length(S)
% 
%     radial_rms_SILVER{n} = zeros(1,64);
%     radial_max_SILVER{n} = zeros(1,64);
%     
%     radial_rms_GR{n} = zeros(1,64);
%     radial_max_GR{n} = zeros(1,64);
%     
%     radial_rms_UNIFORM{n} = zeros(1,64);
%     radial_max_UNIFORM{n} = zeros(1,64);
%     
%     for angle = 0:0.1:179.9
%         psf_rotated = imrotate(SILVER_PSF{n}, angle, 'bilinear', 'crop');
%         line_squared = psf_rotated(round(end/2),:).^2;
%         radial_rms_SILVER{n} = radial_rms_SILVER{n}+abs(line_squared);
%         radial_max_SILVER{n} = max(radial_max_SILVER{n}, abs(sqrt(line_squared)));
%         
%         psf_rotated = imrotate(GR_PSF{n}, angle, 'bilinear', 'crop');
%         line_squared = psf_rotated(round(end/2),:).^2;
%         radial_rms_GR{n} = radial_rms_GR{n}+abs(line_squared);
%         radial_max_GR{n} = max(radial_max_GR{n}, abs(sqrt(line_squared)));
%         
%         psf_rotated = imrotate(Uniform_PSF{n}, angle, 'bilinear', 'crop');
%         line_squared = psf_rotated(round(end/2),:).^2;
%         radial_rms_UNIFORM{n} = radial_rms_UNIFORM{n}+abs(line_squared);
%         radial_max_UNIFORM{n} = max(radial_max_UNIFORM{n}, abs(sqrt(line_squared)));
% 
%     end
% 
%     radial_rms_SILVER{n} = sqrt(radial_rms_SILVER{n}./1800);
%     radial_rms_GR{n} = sqrt(radial_rms_GR{n}./1800);
%     radial_rms_UNIFORM{n} = sqrt(radial_rms_UNIFORM{n}./1800);
%     
%     figure(n+1)
%     subplot(1,2,1)
%     title('RMS')
%     plot(log(radial_rms_SILVER{n}), 'linewidth', 4);
%     set(gca,'FontSize',24)
%     hold on
%     plot(log(radial_rms_GR{n}), 'linewidth', 4);
%     plot(log(radial_rms_UNIFORM{n}), 'linewidth', 4);
% 
%     subplot(1,2,2)
%     title('Max')
%     plot(log(radial_max_SILVER{n}), 'linewidth', 4);
%     set(gca,'FontSize',24)
%     hold on
%     plot(log(radial_max_GR{n}), 'linewidth', 4);
%     plot(log(radial_max_UNIFORM{n}), 'linewidth', 4);
% 
% end

    
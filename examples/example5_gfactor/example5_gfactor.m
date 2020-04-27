%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   How does the trajectory (SILVER,golden ratio, or uniform radial) affect
%   the g-factor maps for a set sensitivity profile?
%                                              
%   Sophie Schauman 2020                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% 1. Choose  set of window sizes, S, to consider
N = [16,32,64,128];
for n = 1:length(N)
    S{n*7-6} = [N(n), 2*N(n)];
    S{n*7-5} = [N(n), 2*N(n), 3*N(n)];
    S{n*7-4} = [N(n)-1:N(n)+1];
    S{n*7-3} = [N(n)-2:N(n)+2];
    S{n*7-2} = [N(n)-3:N(n)+3];
    S{n*7-1} = [N(n)-4:N(n)+4];
    S{n*7} = [N(n)-5:N(n)+5];
end
    

%% 2. Do SILVER optimization for that range

for s = S
   savename = ['examples/precalculated/silver_' strrep(num2str(s{:}),' ', '_') '.mat'];
   if ~exist(savename, 'file')
        SILVER_2D(s{:},'electrostatic_potential',savename) ;
   end
   disp([num2str(s{:}) ' done'])
end

%% 3. Generate g-factors for SILVER

load('examples/example5_gfactor/compressed_sensitivities.mat', 'sens')
psens = zeros(64,64,1,8);
for i = 1:8
    psens(:,:,1,i) = imresize(sens(:,:,1,i),[64,64]);
end

for m = 1:length(S)
    savename = ['examples/example5_gfactor/example5_gfactor_SILVER_' strrep(num2str(S{m}),' ', '_') '_result.mat'];
    if ~exist(savename, 'file')
        for n = 1:length(S{m})
            load(['examples/precalculated/silver_' strrep(num2str(S{m}),' ', '_') '.mat'],'ratio');
            k_SILVER= reshape(sqrt(2)*gen_radial_traj((0:S{m}(n)-1)*ratio*pi, 128, [])./pi,[],2);
            SILVER_g{n} = gFactor(psens,k_SILVER);
        end
        save(savename, 'SILVER_g')
        SILVER_gmaps{m} = SILVER_g;
        clear SILVER_g
    else
        load(savename, 'SILVER_g')
        warning('Using precalculated g-factor maps)')
        SILVER_gmaps{m} = SILVER_g;
        clear SILVER_g
    end
end

%% 4. Generate g-factors for golden ratio and uniform
for n = 1:length(N)
    savename = ['examples/example5_gfactor/example5_gfactor_GR_' num2str(N(n)) '_result.mat'];
    if ~exist(savename, 'file')
        k_GR= reshape(sqrt(2)*gen_radial_traj((0:N(n)-1)*gr2D*pi, 128, [])./pi,[],2);
        GR_g = gFactor(psens,k_GR);
        save(savename, 'GR_g')
        GR_gmaps{n} = GR_g;
        clear GR_g
    else
        load(savename, 'GR_g')
        warning('Using precalculated g-factor maps)')
        GR_gmaps{n} = GR_g;
        clear GR_g
    end
    savename = ['examples/example5_gfactor/example5_gfactor_UNIFORM_' num2str(N(n)) '_result.mat'];
    if ~exist(savename, 'file')
        k_UNIFORM= reshape(sqrt(2)*gen_radial_traj((0:N(n)-1)*1/N(n)*pi, 128, [])./pi,[],2);
        UNIFORM_g = gFactor(psens,k_UNIFORM);
        save(savename, 'UNIFORM_g')
        UNIFORM_gmaps{n} = UNIFORM_g;
        clear UNIFORM_g
    else
        load(savename, 'UNIFORM_g')
        warning('Using precalculated g-factor maps)')
        UNIFORM_gmaps{n} = UNIFORM_g;
        clear UNIFORM_g
    end
    
end


%% 5. Visualize
[columnsInImage, rowsInImage] = meshgrid(1:64, 1:64);
centerX = 32;
centerY = 32;
radius = 30;
circleMask = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;

p = 1;
v = repmat(1:length(N),[7 1]);
v = v(:);
for m = 1:length(S)
    a = ismember(S{m},N(v(m)));
    S{m}(a)
    subplot(length(N),7,p)
    im = abs(SILVER_gmaps{m}{a}).*circleMask;
    SILVER_max_g(m)= max(im(circleMask));
    SILVER_mean_g(m)= mean(im(circleMask));
    imagesc(im)
    axis image
    axis off
    if p == 1
        title('S = \{N, 2N\}')
    elseif p == 2
        title('S = \{N, 2N, 3N\}')
    elseif p == 3
        title('S = \{N-1 ... N+1\}')
    elseif p == 4
        title('S = \{N-2 ... N+2\}')
    elseif p == 5
        title('S = \{N-3 ... N+3\}')
    elseif p == 6
        title('S = \{N-4 ... N+4\}')
    elseif p == 7
        title('S = \{N-5 ... N+5\}')
    end
    
    if mod(p,7) == 1 
        cl = [min(im(im>0)),max(im(:))];
        caxis(cl);
        text(-2,32, ['N = ' num2str(N(v(m))) ' spokes'], 'fontsize', 14,'HorizontalAlignment','right')
    else
        caxis(cl)
    end
    p = p+1;
end

set(gcf, 'Position', [440 86 888 712])

SILVER_max_g= reshape(SILVER_max_g,7,length(N));
SILVER_mean_g= reshape(SILVER_mean_g,7,length(N));

%% 6. Compare SILVER to GR and Uniform
lables = {'S = \{N, 2N\}', 'S = \{N, 2N, 3N\}', 'S = \{N-1 ... N+1\}', 'S = \{N-2 ... N+2\}','S = \{N-3 ... N+3\}','S = \{N-4 ... N+4\}','S = \{N-5 ... N+5\}', 'UNIFORM', 'Golden ratio'};
for m = 1:length(N)
    figure
    hold on
    [data, idx] = sort(abs(cat(1,SILVER_max_g(:,m),max(UNIFORM_gmaps{m}(circleMask)),max(GR_gmaps{m}(circleMask)))));
    for ii = 1:length(data)
        h = bar(ii,data(ii));
        if idx(ii) == 8 %UNIFORM
            set(h, 'FaceColor', 'r') 
        elseif idx(ii) == 9 %golden ratio
            set(h, 'FaceColor', 'y') 
        else %SILVER
            set(h, 'FaceColor', 'b') 
        end
    end
    xticks(1:9)
    xticklabels(lables(idx))
    xtickangle(-90)
    title(['N = ' num2str(N(m))])
end





% %% 4. Visualize
% figure(1)
% [columnsInImage, rowsInImage] = meshgrid(1:64, 1:64);
% centerX = 32;
% centerY = 32;
% radius = 60;
% circleMask = (rowsInImage - centerY).^2 ...
%     + (columnsInImage - centerX).^2 <= radius.^2;
% 
% c = [1,100; 1,10; 1, 1.1];
% 
% for n = 1:length(S)
%     subplot(length(S),3,n*3-2)
%     im = abs(Uniform_g{n}).*circleMask;
%     disp(['Maximum g-factor for Uniform sampling (' num2str(S(n)) ' spokes) is ' num2str(max(im(:)))])
%     imagesc(im)
%     axis image
%     axis off
%     if n ==1
%         text(32,-4,'UNIFORM', 'fontsize', 20,'HorizontalAlignment','center')
%     end
%     caxis(c(n,:))
%     colorbar
%     set(gca,'FontSize', 18)
%     text(-2,32, [num2str(S(n)) ' spokes'], 'fontsize', 20,'HorizontalAlignment','right')
%     
%     
%     subplot(length(S),3,n*3-1)
%     im = abs(GR_g{n}).*circleMask;
%     disp(['Maximum g-factor for Golden ratio sampling (' num2str(S(n)) ' spokes) is ' num2str(max(im(:)))])
%     imagesc(im)
%     axis image
%     axis off
%     if n ==1
%         text(32,-4,'GOLDEN RATIO', 'fontsize', 20,'HorizontalAlignment','center')
%     end
%     caxis(c(n,:))
%     colorbar
%     set(gca,'FontSize', 18)
%     
%     subplot(length(S),3,n*3)
%     im = abs(SILVER_g{n}).*circleMask;
%     disp(['Maximum g-factor for SILVER sampling (' num2str(S(n)) ' spokes) is ' num2str(max(im(:)))])
%     imagesc(im)
%     axis image
%     axis off
%     if n ==1
%         text(32,-4,'SILVER', 'fontsize', 20,'HorizontalAlignment','center')
%     end
%     caxis(c(n,:));
%     colorbar
%     set(gca,'FontSize', 18)
%     disp(' ')
%     
% end
% % set(gcf,'Position', [440 1 794 797])
% % savefig('examples/example5_gfactor/example5_gfactor_result.fig')
% % saveas(gcf,'examples/example5_gfactor/example5_gfactor_result.tiff')

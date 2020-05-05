%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   How does the trajectory (SILVER,golden ratio, or uniform radial) affect
%   the g-factor maps for a set sensitivity profile?
%                                              
%   Sophie Schauman 2020                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

savefolder = 'examples/example5_gfactor/gfactor_maps/test/';
load('examples/example5_gfactor/test_sensitivity_map.mat', 'psens') %coils

%% 1. Choose set of window sizes, S, to consider
N = [16,32,48,64,128];
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

%% 3. Generate g-factor maps for GOLDEN RATIO and UNIFORM to use for comparison

for m = 1:length(S)
    savename = [savefolder 'example5_gfactor_SILVER_' strrep(num2str(S{m}),' ', '_') '_result.mat'];
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
    savename = [savefolder 'example5_gfactor_GR_' num2str(N(n)) '_result.mat'];
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
    savename = [savefolder 'example5_gfactor_UNIFORM_' num2str(N(n)) '_result.mat'];
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


%% 5. Visual check and preprocess for analysis

figure(1)
p = 1;
v = repmat(1:length(N),[7 1]);
v = v(:);
for m = 1:length(S)
    a = ismember(S{m},N(v(m)));
    b = ismember(N,N(v(m)));
    S{m}(a);
    N(b);
    subplot(length(N),7,p)
    im = cat(2, abs(UNIFORM_gmaps{b}),abs(SILVER_gmaps{m}{a}), abs(GR_gmaps{b}));
    SILVER_max_g(m)= max(abs(SILVER_gmaps{m}{a}(:)));
    SILVER_mean_g(m)= mean(abs(SILVER_gmaps{m}{a}(:)));
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
    data(m,:) = abs(cat(1,SILVER_mean_g(:,m),mean(UNIFORM_gmaps{m}(:)),mean(GR_gmaps{m}(:))));
    data(m,:) = data(m,:)./abs(mean(UNIFORM_gmaps{m}(:)));
end
  
figure(2)
hold on

data = mean(data,3); % across pixels

data_means = mean(data,1); % across different N
data_min = data_means-min(data,[],[1,3]);
data_max = max(data,[],[1,3])-data_means;

[data_means, idx] = sort(data_means);
data_min = data_min(idx);
data_max = data_max(idx);
data = data(:,idx);

silver_cl = linspace(0.4,0.7,7)'.*ones(7,3);
for ii = 1:length(data_means)
    h = bar(ii,data_means(ii));
    if idx(ii) == 8 %UNIFORM
        set(h, 'FaceColor', [0 0.5 1]) 
    elseif idx(ii) == 9 %golden ratio
        set(h, 'FaceColor', [1 0.5 0]) 
    else %SILVER
        set(h, 'FaceColor', silver_cl(idx(ii),:)) 
    end
    h.LineWidth = 2;
end
% er = errorbar(1:9,data_means,data_min,data_max);
% er.Color = [0 0 0];                            
% er.LineStyle = 'none'; 
% er.LineWidth = 2;
cl = hsv(length(N));
for n = 1:length(N)
    h(n) = scatter(1:9, data(n,:),100, cl(n,:),'filled','MarkerFaceAlpha', 0.5, 'markeredgecolor','k', 'DisplayName',['N = ' num2str(N(n))]);
end
l = legend(h);
l.Location = 'northwest';
grid on
xticks(1:9)
xticklabels(lables(idx))
xtickangle(-90)
ylabel('mean g-factor compared to uniform')
set(gca,'fontsize', 16)
set(gca,'linewidth', 2)
box on
savefig([savefolder 'example5_gfactor_graph.fig'])
saveas(gcf,[savefolder '/example5_gfactor_graph.tiff'])

%% 7. Show specific example
% choose sets:
%    1 = {N, 2N}
%    2 = {N, 2N, 3N}
%    3 = {N-1, ..., N+1}
%    4 = {N-2, ..., N+2}
%    5 = {N-3, ..., N+3}
%    6 = {N-4, ..., N+4}
%    7 = {N-5, ..., N+5}

chosen_set = 1:7; % chosen set
for cs = chosen_set
    sets = [];
    for n = 1:length(N)
        sets = cat(1,sets,n*7-(7-cs));
    end

    ii = 1;
    im = [];
    for s = S(sets)
        a = ismember(s{:},N(ii));
        s{:}(a)
        im = cat(2,im,cat(1, abs(SILVER_gmaps{sets(ii)}{a})./abs(UNIFORM_gmaps{ii}), abs(GR_gmaps{ii})./abs(UNIFORM_gmaps{ii})));
        ii = ii+1;
    end
    figure
    imagesc(im)
    cl = caxis;
    caxis([1 cl(2)*0.5])
    colormap('jet')
    colorbar
    axis image
    axis off
    text(-max(size(psens)),max(size(psens))/2,{'SILVER:' lables{cs}}, 'fontsize', 16,'HorizontalAlignment','left')
    text(-max(size(psens)),max(size(psens))+max(size(psens))/2,'Golden ratio', 'fontsize', 16,'HorizontalAlignment','left')
    for n = 1:length(N)
        text((n-1)*max(size(psens))+max(size(psens))/2,-max(size(psens))/8,['N = ' num2str(N(n))], 'fontsize', 16,'HorizontalAlignment','center')
    end
    set(gcf, 'Position', [57 444 1167 354])

    savefig([savefolder 'example5_gfactor_exampleset' num2str(cs) '.fig'])
    saveas(gcf,[savefolder 'example5_gfactor_exampleset' num2str(cs) '.tiff'])
end
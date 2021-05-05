function [] = measure_noise_amplification_multi_coil_fig(spokes,sensmask,data_files, savename)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


a = dir(data_files);

n=0;
for ii = 1:length(a)
    if contains(a(ii).name, '_GR_') 
        load(a(ii).name, 'GR_n');
    elseif contains(a(ii).name, '_U_')
        load(a(ii).name, 'U_n');
    elseif contains(a(ii).name, '_SILVER_')        
        n = n+1;
        nums = regexp(a(ii).name,'\d*','Match');
        S_set{n} = [str2double(nums(3:end))];
        S_rats(n) = SILVER_2D(S_set{n}, 'electrostatic_potential');
        tmp = load(a(ii).name, 'SILVER_n');
        S_n(:,:,n) = tmp.SILVER_n;
        clear tmp
    end
end
 
col =[1 0.5 0; 0 0.5 1; ones(size(S_n,3),3).*linspace(0.4, 0.7, size(S_n,3))'];

subplot(3,5,1)
show_spokes(gr2D, [0:spokes-1])
axis image
xticks([])
yticks([])
set(gca, 'XColor', col(1,:))
set(gca, 'YColor', col(1,:))
set(gca, 'LineWidth', 2)
title('GR', 'FontSize', 12)
text(-1.5, 1.5, '(A)', 'fontsize', 14);



subplot(3,5,2)
show_spokes(1/spokes, [0:spokes-1])
axis image
xticks([])
yticks([])
set(gca, 'XColor', col(2,:))
set(gca, 'YColor', col(2,:))
set(gca, 'LineWidth', 2)
title('Uniform', 'FontSize', 12)

subplot(3,5,3)
show_spokes(S_rats(1), [0:spokes-1])
axis image
xticks([])
yticks([])
set(gca, 'XColor', col(3,:))
set(gca, 'YColor', col(3,:))
set(gca, 'LineWidth', 2)
title(['SILVER \{' num2str(S_set{1}) '\}'], 'FontSize', 12)

subplot(3,5,4)
show_spokes(S_rats(2), [0:spokes-1])
axis image
xticks([])
yticks([])
set(gca, 'XColor', col(4,:))
set(gca, 'YColor', col(4,:))
set(gca, 'LineWidth', 2)
title(['SILVER \{' num2str(S_set{2}) '\}'],'FontSize', 12)

subplot(3,5,5)
show_spokes(S_rats(3), [0:spokes-1])
axis image
xticks([])
yticks([])
set(gca, 'XColor', col(5,:))
set(gca, 'YColor', col(5,:))
set(gca, 'LineWidth', 2)
title(['SILVER \{' num2str(S_set{3}) '\}'], 'FontSize', 12)


subplot(3,5,[6:10])
imagesc(cat(2, GR_n.*sensmask, U_n.*sensmask,reshape(S_n, size(S_n,1),[]).*repmat(sensmask,1,length(S_rats))))
axis image
axis off
text(15, -3, 'GR', 'HorizontalAlignment', 'center', 'FontSize', 12)
text(15+30, -3, 'Uniform', 'HorizontalAlignment', 'center', 'FontSize', 12)
text(15+30*2, -3, ['SILVER \{' num2str(S_set{1}) '\}'], 'HorizontalAlignment', 'center', 'FontSize', 12)
text(15+30*3, -3, ['SILVER \{' num2str(S_set{2}) '\}'], 'HorizontalAlignment', 'center', 'FontSize', 12)
text(15+30*4, -3, ['SILVER \{' num2str(S_set{3}) '\}'], 'HorizontalAlignment', 'center', 'FontSize', 12)
text(-3, 15, 'noise s.d. maps', 'HorizontalAlignment', 'center', 'FontSize', 12, 'Rotation', 90)

text(-5, -3, '(B)', 'fontsize', 14);


colormap hot
colorbar
subplot(3,2,5)
[U_GR_p, U_GR_h] = signrank(GR_n(sensmask), U_n(sensmask), 'alpha', 0.05/(1+2*size(S_n,3)+nchoosek(size(S_n,3),2)));

mean_noise = [mean(GR_n(sensmask), 'all'), mean(U_n(sensmask), 'all')];
std_noise = [std(GR_n(sensmask), [],'all'), std(U_n(sensmask), [],'all')];
eff = [efficiency_2D(gr2D*pi*[0:spokes-1], 'electrostatic_potential'),1];
for ii=1:size(S_n,3)
    tmp = S_n(:,:,ii);
    mean_noise(ii+2) = mean(tmp(sensmask), 'all');
    std_noise(ii+2) = std(tmp(sensmask), [],'all');
    eff(ii+2) = efficiency_2D(S_rats(ii)*pi*[0:spokes-1], 'electrostatic_potential');
    [U_S_p(ii), U_S_h(ii)] = signrank(U_n(sensmask), tmp(sensmask), 'alpha', 0.05/(1+2*size(S_n,3)+nchoosek(size(S_n,3),2)));
    [GR_S_p(ii), GR_S_h(ii)] = signrank(GR_n(sensmask), tmp(sensmask), 'alpha', 0.05/(1+2*size(S_n,3)+nchoosek(size(S_n,3),2)));
end

S_S_h = nan(size(S_n,3)-1,size(S_n,3));
for n = 1:size(S_n,3)-1
    for m = n+1:size(S_n,3)
        tmp_n = S_n(:,:,n);
        tmp_m = S_n(:,:,m);
        [S_S_p(n,m), S_S_h(n,m)] = signrank(tmp_n(sensmask), tmp_m(sensmask), 'alpha', 0.05/(1+2*size(S_n,3)+nchoosek(size(S_n,3),2)));
    end
end

for n = 1:size(col,1)
    b = bar(n,mean_noise(n));
    hold on
    b.FaceColor = 'Flat';
    b.CData = col(n,:);
end
errorbar(1:length(mean_noise),mean_noise,std_noise, 'LineStyle', 'none', 'color', 'k')


if ~U_GR_h
    line([1 2], [max(mean_noise(1:2)+std_noise(1:2))*1.03, max(mean_noise(1:2)+std_noise(1:2))*1.03], 'color', 'k', 'handlevisibility', 'off', 'marker', '.')
    text([1.5], [max(mean_noise(1:2)+std_noise(1:2))*1.04], 'n.s.', 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
end
for n = 1:size(S_n,3)
    if ~GR_S_h(n)
        line([1 n+2], [max(mean_noise([1,n+2])+std_noise([1,n+2]))*1.02, max(mean_noise([1,n+2])+std_noise([1,n+2]))*1.02], 'color', 'k', 'handlevisibility', 'off', 'marker', '.')
        text([(1+n+2)/2], [max(mean_noise([1,n+2])+std_noise([1,n+2]))*1.03], 'n.s.', 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
    end
    if ~U_S_h(n)
        line([2 n+2], [max(mean_noise([2,n+2])+std_noise([2,n+2]))*1.02, max(mean_noise([2,n+2])+std_noise([2,n+2]))*1.02], 'color', 'k', 'handlevisibility', 'off', 'marker', '.')
        text([(2+n+2)/2], [max(mean_noise([2,n+2])+std_noise([2,n+2]))*1.03], 'n.s.', 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
    end
    for m = n+1:size(S_n,3)
        if ~S_S_h(n,m)
            line([n+2 m+2], [max(mean_noise([m+2,n+2])+std_noise([m+2,n+2]))*1.02, max(mean_noise([m+2,n+2])+std_noise([m+2,n+2]))*1.02], 'color', 'k', 'handlevisibility', 'off', 'marker', '.')
            text([(n+2+m+2)/2], [max(mean_noise([m+2,n+2])+std_noise([m+2,n+2]))*1.03], 'n.s.', 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        end
    end
        
end


axis([xlim min(mean_noise(:)-std_noise(:))*0.99  max(mean_noise(:)+std_noise(:))*1.05])
ylabel('noise s.d.', 'Color', 'k')
legend('GR', 'Uniform', ['SILVER \{' num2str(S_set{1}) '\}'],['SILVER \{' num2str(S_set{2}) '\}'],['SILVER \{' num2str(S_set{3}) '\}'],'Location', 'southoutside')
xticks([])
grid minor
text(-1, max(ylim), '(C)', 'FontSize',14, 'VerticalAlignment', 'bottom')

subplot(3,2,6)
for n = 1:size(col,1)
    b = bar(n,1./eff(n));
    hold on
    b.FaceColor = 'Flat';
    b.CData = col(n,:);
end

axis([xlim min(1./eff(:))*0.99  max(1./eff(:))*1.01])
ylabel('1/efficiency', 'Color', 'k')
legend('GR', 'Uniform', ['SILVER \{' num2str(S_set{1}) '\}'],['SILVER \{' num2str(S_set{2}) '\}'],['SILVER \{' num2str(S_set{3}) '\}'],'Location', 'southoutside')
xticks([])
grid minor
text(-1, max(ylim), '(D)', 'FontSize',14, 'VerticalAlignment', 'bottom')

annotation('line', [0 1], [0.7 0.7]);
annotation('line', [0 1], [0.37 0.37]);
annotation('line', [0.5 0.5], [0.37 0]);
set(gcf, 'Position', [141 71 1019 734])
savefig(gcf, [savename '.fig'])
saveas(gcf, [savename '.svg'])

end


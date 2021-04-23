function [] = measure_noise_amplification_single_coil_fig(spokes,data_files, savename)
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
        S_set{n} = [str2double(nums(2:end))];
        S_rats(n) = SILVER_2D(S_set{n}, 'electrostatic_potential');
        tmp = load(a(ii).name, 'SILVER_n');
        S_n(:,:,n) = tmp.SILVER_n;
        clear tmp
    end
end
 

col =[1 0.5 0; 0 0.5 1; ones(size(S_n,3),3).*linspace(0.4, 0.7, size(S_n,3))'];

subplot(2,5,1)
show_spokes(gr2D, [0:spokes-1])
axis image
xticks([])
yticks([])
set(gca, 'XColor', col(1,:))
set(gca, 'YColor', col(1,:))
set(gca, 'LineWidth', 2)
title('GR', 'FontSize', 12)
text(-1.5, 1.5, '(A)', 'fontsize', 14);

subplot(2,5,2)
show_spokes(1/spokes, [0:spokes-1])
axis image
xticks([])
yticks([])
set(gca, 'XColor', col(2,:))
set(gca, 'YColor', col(2,:))
set(gca, 'LineWidth', 2)
title('Uniform', 'FontSize', 12)

subplot(2,5,3)
show_spokes(S_rats(1), [0:spokes-1])
axis image
xticks([])
yticks([])
set(gca, 'XColor', col(3,:))
set(gca, 'YColor', col(3,:))
set(gca, 'LineWidth', 2)
title(['SILVER \{' num2str(S_set{1}) '\}'], 'FontSize', 12)

subplot(2,5,4)
show_spokes(S_rats(2), [0:spokes-1])
axis image
xticks([])
yticks([])
set(gca, 'XColor', col(4,:))
set(gca, 'YColor', col(4,:))
set(gca, 'LineWidth', 2)
title(['SILVER \{' num2str(S_set{2}) '\}'],'FontSize', 12)

subplot(2,5,5)
show_spokes(S_rats(3), [0:spokes-1])
axis image
xticks([])
yticks([])
set(gca, 'XColor', col(5,:))
set(gca, 'YColor', col(5,:))
set(gca, 'LineWidth', 2)
title(['SILVER \{' num2str(S_set{3}) '\}'], 'FontSize', 12)


subplot(2,2,3)  
mean_noise = [mean(GR_n, 'all'), mean(U_n, 'all')];
std_noise = [std(GR_n, [],'all'), std(U_n, [],'all')];
eff = [efficiency_2D(gr2D*pi*[0:spokes-1], 'electrostatic_potential'),1];
for ii=1:size(S_n,3)
    mean_noise(ii+2) = mean(S_n(:,:,ii), 'all');
    std_noise(ii+2) = std(S_n(:,:,ii), [],'all');
    eff(ii+2) = efficiency_2D(S_rats(ii)*pi*[0:spokes-1], 'electrostatic_potential');
end


for n = 1:size(col,1)
    b = bar(n,mean_noise(n));
    hold on
    b.FaceColor = 'Flat';
    b.CData = col(n,:);
end
errorbar(1:length(mean_noise),mean_noise,std_noise, 'LineStyle', 'none', 'color', 'k')

axis([xlim min(mean_noise(:)-std_noise(:))*0.99  max(mean_noise(:)+std_noise(:))*1.01])
ylabel('noise', 'Color', 'k')
legend('GR', 'Uniform', ['SILVER \{' num2str(S_set{1}) '\}'],['SILVER \{' num2str(S_set{2}) '\}'],['SILVER \{' num2str(S_set{3}) '\}'],'Location', 'southoutside')
xticks([])
grid minor
text(-1, max(ylim), '(B)', 'FontSize',14, 'VerticalAlignment', 'bottom')


subplot(2,2,4)
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
text(-1, max(ylim), '(C)', 'FontSize',14, 'VerticalAlignment', 'bottom')

annotation('line', [0 1], [0.55 0.55]);

set(gcf, 'Position', [141 164 998 495])
savefig(gcf, [savename '.fig'])
saveas(gcf, [savename '.svg'])

end


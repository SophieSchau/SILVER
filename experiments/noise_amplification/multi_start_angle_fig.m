function [] = multi_start_angle_fig(start_angles,sensmask, data_files, savename)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


for start_angle = [0, start_angles]
    a = dir( [data_files num2str(start_angle) '_*']);

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
            tmp = load(a(ii).name, 'SILVER_n');
            S_n(:,:,n) = tmp.SILVER_n;
            clear tmp
        end
    end
    GR(start_angle+1) = mean(GR_n(sensmask),'all');
    U(start_angle+1) = mean(U_n(sensmask),'all');
    for nn = 1:size(S_n,3)
        tmp = S_n(:,:,nn);
        S(start_angle+1, nn) = mean(tmp(sensmask),'all');
    end

end

% significance testing
U_GR_ttest = ttest(GR(:), U(:), 'alpha', 0.05/(1+2*size(S_n,3)+nchoosek(size(S_n,3),2)));
for ii=1:size(S_n,3)
    U_S_ttest(ii) = ttest(U(:), S(:,ii), 'alpha', 0.05/(1+2*size(S_n,3)+nchoosek(size(S_n,3),2)));
    GR_S_ttest(ii) = ttest(GR(:), S(:,ii), 'alpha', 0.05/(1+2*size(S_n,3)+nchoosek(size(S_n,3),2)));
end
S_S_ttest = nan(size(S_n,3)-1,size(S_n,3));
for n = 1:size(S_n,3)-1
    for m = n+1:size(S_n,3)
        tmp_n = S(:,n);
        tmp_m = S(:,m);
        S_S_ttest(n,m) = ttest(tmp_n(:), tmp_m(:), 'alpha', 0.05/(1+2*size(S_n,3)+nchoosek(size(S_n,3),2)));
    end
end


col =[1 0.5 0; 0 0.5 1; ones(size(S_n,3),3).*linspace(0.4, 0.7, size(S_n,3))'];

plot([0,start_angles], GR, 'color',col(1,:), 'linewidth', 2);
hold on

plot([0,max(start_angles)], sqrt([mean(GR.^2) mean(GR.^2)]), 'color',col(1,:), 'linewidth', 2, 'LineStyle', ':', 'HandleVisibility', 'off');

plot([0,start_angles], U, 'color',col(2,:),'linewidth', 2);
plot([0,max(start_angles)], sqrt([mean(U.^2) mean(U.^2)]), 'color',col(2,:), 'linewidth', 2, 'LineStyle', ':', 'HandleVisibility', 'off');

for ss = 1:size(S_n,3)
    plot([0,start_angles], S(:,ss), 'color',col(ss+2,:),'linewidth', 2);
    plot([0,max(start_angles)], sqrt([mean(S(:,ss).^2) mean(S(:,ss).^2)]), 'color',col(ss+2,:), 'linewidth', 2, 'LineStyle', ':', 'HandleVisibility', 'off');
end

if ~U_GR_ttest
    line([start_angles(end)+1 start_angles(end)+1], [sqrt(mean(GR.^2)) sqrt(mean(U.^2))], 'color', 'k', 'handlevisibility', 'off', 'marker', '.')
    text([start_angles(end)+1], mean([sqrt(mean(GR.^2)) sqrt(mean(U.^2))]), 'n.s.', 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', -90)
end
for n = 1:size(S_n,3)
    if ~GR_S_ttest(n)
        line([start_angles(end)+1 start_angles(end)+1], [sqrt(mean(GR.^2)) sqrt(mean(S(:,n).^2))], 'color', 'k', 'handlevisibility', 'off', 'marker', '.')
        text([start_angles(end)+1], mean([sqrt(mean(GR.^2)) sqrt(mean(S(:,n).^2))]), 'n.s.', 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', -90)
    end
    if ~U_S_ttest(n)
        line([start_angles(end)+1 start_angles(end)+1], [sqrt(mean(U.^2)) sqrt(mean(S(:,n).^2))], 'color', 'k', 'handlevisibility', 'off', 'marker', '.')
        text([start_angles(end)+1], mean([sqrt(mean(U.^2)) sqrt(mean(S(:,n).^2))]), 'n.s.', 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', -90)
    end
    for m = n+1:size(S_n,3)
        if ~S_S_ttest(n,m)
            line([start_angles(end)+1 start_angles(end)+1], [sqrt(mean(S(:,m).^2)) sqrt(mean(S(:,n).^2))], 'color', 'k', 'handlevisibility', 'off', 'marker', '.')
            text([start_angles(end)+1], mean([sqrt(mean(S(:,m).^2)) sqrt(mean(S(:,n).^2))]), 'n.s.', 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', -90)
        end
    end
end



xlabel('start angle, degrees')
ylabel('noise')
set(gca, 'fontsize', 14)
legend('GR', 'Uniform', ['SILVER \{' num2str(S_set{1}) '\}'],['SILVER \{' num2str(S_set{2}) '\}'],['SILVER \{' num2str(S_set{3}) '\}'],'Location', 'southoutside')
set(gcf, 'Position', [436 265 881 413])
grid minor
box off
savefig(gcf, [savename '.fig'])
saveas(gcf, [savename '.svg'])



end


function [] = make_nifti_from_slices(file_list,resolutions,slices,savename)

    for ii = 1:length(file_list)
        load(file_list(ii).name);
        sl_idx1 = regexp(file_list(ii).name, '_sl._');
        sl_idx2 = regexp(file_list(ii).name, '_sl.._');

        if ~isempty(sl_idx1)
            sl = str2double(file_list(ii).name(sl_idx1+3));
        elseif ~isempty(sl_idx2)
            sl = str2double(file_list(ii).name(sl_idx2+3:sl_idx2+4));
        end

        recon_vol(:,:,sl,:) = recon_a;
    end
    recon_vol = imrotate(recon_vol(:,:,slices,:), 180);
    save_avw(abs(recon_vol),savename,'d',resolutions)
end
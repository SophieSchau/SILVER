function [] = read_in_and_save_TURBINE(raw_filepath,savepath)
%READ_IN_AND_SAVE_TURBINE Prepares TURBINE data for 2D reconstruction
%   This function takes all .dat TURBINE files in a folder, performs nyquist
%   ghost correction as well as global phase correction. It assumes that the
%   data is fully sampled in the PE direction (z) and performs an iFFT along
%   z. It then saves the slices (without metadata) into separate folders.
%
%   Sophie Schauman 2020
%%%%%%%%%%%%%

    files = dir(raw_filepath);
    for scan = [1:length(files)]
        raw_filename = files(scan).name;
        
        if contains(raw_filename, '.dat')

            twix_obj = mapVBVD([raw_filepath, raw_filename], 'removeOS', false,'rampsampregrid',true);

            [kdata_3D,pcorr_nyquist] = phaseCorr(twix_obj{2}.phasecor(), twix_obj{2}.image(), 'mode', 'TURBINE');
            [kdata_3D,pcorr_shot] = physioCorr(reshape(kdata_3D,size(kdata_3D,1), size(kdata_3D,2), size(kdata_3D,3),[]));

            % Fourier transform accross z (fully sampled)
            kdata_2D = ifftdim(kdata_3D,3);

            for slice = 1:size(kdata_3D,3)
                kdata = permute(kdata_2D(:,:,slice,:),[1,4,3,2]);
                save([savepath raw_filename(1:end-4) '_slice' num2str(slice)], 'kdata', 'pcorr_nyquist', 'pcorr_shot')
            end
        end
    end
end


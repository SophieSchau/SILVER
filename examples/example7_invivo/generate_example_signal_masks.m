%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The signal masks to use for quality assessment.
%
%   Sophie Schauman 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

load('examples/example7_invivo/example_kdata_68_153_306.mat', 'kdata_Uniform', 'S', 'SILVER_twix');
load('examples/example7_invivo/example_sensitivities_68_153_306.mat', 'sens_UNIFORM');
[S_max, S_max_idx] = max(S);
NCoils = SILVER_twix.hdr.Meas.NChaMeas;
Mat_size = [SILVER_twix.hdr.Config.BaseResolution, SILVER_twix.hdr.Config.BaseResolution, 1, 1];
%%

ratio = 1/S_max;
NFrames = SILVER_twix.hdr.Meas.NPhs;
NRepeats = SILVER_twix.hdr.Config.NLin/SILVER_twix.hdr.Config.NSeg;
NSpokes = SILVER_twix.hdr.Config.NLin;

Phi = mod([0:(NSpokes*NFrames-1)] * ratio * pi, 2*pi );
kspace = gen_radial_traj(Phi, SILVER_twix.hdr.Config.NColMeas, []);

%%
for subj = 1:size(kdata_Uniform,2)
    E = xfm_NUFFT(Mat_size,sens_UNIFORM{subj},[],...
            reshape(kspace, [], 1, 2));
    kdata_Uniform_presub = reshape(kdata_Uniform{S_max_idx,subj}(:,:,2,:)-kdata_Uniform{S_max_idx,subj}(:,:,1,:),[],1,NCoils);

    recon = fista(E, kdata_Uniform_presub, 1, 0.0000001, ...
        Mat_size, 10, 0.5);

    
    mask_signal{subj} = imbinarize(abs(recon),0.0002);
    mask_noise{subj} = logical(logical(abs(recon)).*~imbinarize(abs(recon),0.0001));
    
    figure
    imagesc(abs(recon))
    hold on
    h = imagesc(mask_signal{subj});
    set(h,'AlphaData', 0.5)
    h2 = imagesc(mask_noise{subj});
    set(h2,'AlphaData', 0.5)
    drawnow

end
    
save('examples/example7_invivo/example_masks_68_153_306.mat', 'mask_signal', 'mask_noise', '-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The signal masks to use for quality assessment.
%
%   Sophie Schauman 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
load('examples/example7_invivo/example_kdata_UNIFORM_68_153_306.mat','kdata_Uniform');
load('examples/example7_invivo/example_params_68_153_306.mat');
load('examples/example7_invivo/example_sensitivities_68_153_306.mat', 'sens_UNIFORM');
[S_max, S_max_idx] = max(S);
%%

ratio = 1/S_max;

Phi = mod([0:(NSpokes*NFrames-1)] * ratio * pi, 2*pi );
kspace = gen_radial_traj(Phi, NSamps, []);

%%
for subj = 1:size(kdata_Uniform,2)
    E = xfm_NUFFT([Mat_size, Mat_size, 1, 1],sens_UNIFORM{subj},[],...
            reshape(kspace, [], 1, 2));
    kdata_Uniform_presub = reshape(kdata_Uniform{S_max_idx,subj}(:,:,2,:)-kdata_Uniform{S_max_idx,subj}(:,:,1,:),[],1,size(sens_UNIFORM{1,subj},4));

    recon = fista(E, kdata_Uniform_presub, 1, 0.0000001, ...
        [Mat_size, Mat_size, 1, 1], 10, 0.5);

    BW = bwareafilt(imclose(imbinarize(abs(recon*100)),ones(4)),1);
    s = regionprops(BW,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
    X0=s.Centroid(1); %Coordinate X
    Y0=s.Centroid(2); %Coordinate Y
    l=s.MajorAxisLength*0.9/2; %Length
    w=s.MinorAxisLength*0.9/2; %Width
    phi=s.Orientation; 
    [X, Y] = meshgrid(1:size(recon,1),1:size(recon,2)); 
    ellipse = ((X-X0)/w).^2+((Y-Y0)/l).^2<=1;
    mask = ellipse;
    mask_signal{subj} = imbinarize(abs(recon),0.0002);
    mask_noise{subj} = logical(mask.*~imbinarize(abs(recon),0.0001));

    if subj == 3
        mask_signal{subj} = imbinarize(abs(recon),0.0003);
        mask_noise{subj} = logical(mask.*~imbinarize(abs(recon),0.0002));

    end
    
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
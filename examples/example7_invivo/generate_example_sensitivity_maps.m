%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The example maps here were generated by time-averaging ASL Angiography
%   data generated with 'generate_example_kdata.m'.
%
%   Sophie Schauman 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

load('examples/example7_invivo/example_kdata_68_153_306.mat', 'kdata_SILVER', 'kdata_GR', 'kdata_Uniform', 'S', 'SILVER_twix');
[S_max, S_max_idx] = max(S);
NCoils = SILVER_twix.hdr.Meas.NChaMeas;
Mat_size = [SILVER_twix.hdr.Config.BaseResolution, SILVER_twix.hdr.Config.BaseResolution, 1, 1];
%%

for method = {'UNIFORM', 'GR', 'SILVER'}
    switch method{:}
        case 'UNIFORM'
            ratio = 1/S_max;
        case 'GR'
            ratio = gr2D;
        case 'SILVER'
            ratio = SILVER_2D(S,'electrostatic_potential') ;
    end
    Phi = [];
    for frame = 1:SILVER_twix.hdr.Meas.NPhs
        for repeat = 0:SILVER_twix.hdr.Config.NLin/SILVER_twix.hdr.Config.NSeg-1
            Phi = cat(1, Phi, mod( ((frame-1)*SILVER_twix.hdr.Config.NLin+repeat:SILVER_twix.hdr.Config.NLin/SILVER_twix.hdr.Config.NSeg:frame*SILVER_twix.hdr.Config.NLin-1)' * ratio * pi, 2*pi ));
        end
    end
    switch method{:}
        case 'UNIFORM'
            kspace_UNIFORM = gen_radial_traj(Phi, SILVER_twix.hdr.Config.NColMeas, []);
        case 'GR'
            kspace_GR = gen_radial_traj(Phi, SILVER_twix.hdr.Config.NColMeas, []);
        case 'SILVER'
            kspace_SILVER = gen_radial_traj(Phi, SILVER_twix.hdr.Config.NColMeas, []);
    end
end

%%
for subj = 1:size(kdata_SILVER,2)
    E_UNIFORM = xfm_NUFFT(Mat_size,ones(Mat_size),[],...
            reshape(kspace_UNIFORM, [], 1, 1, 2));
    kdata_Uniform_c = reshape(mean(kdata_Uniform{S_max_idx,subj},3),[],1,1,NCoils);
    
    
    for c = 1:NCoils
        disp(['Calculating coil sensitivity map ' num2str(c) ' of ' num2str(NCoils)])
        ims(:,:,:,c) =   fista(E_UNIFORM, kdata_Uniform_c(:,:,:,c), 1, 0.00000001, ...
            [E_UNIFORM.Nd, 1,1], 10, 0.5);
        imshow(abs(ims(:,:,:,c)),[])
        drawnow
    end
    disp(['Adaptive combine estimation of coil sensitivities  - Uniform sampling'])
    sens_UNIFORM{subj} = adaptive_estimate_sens('data', permute(ims,[4,1,2,3]),...
            'kernel', 10, 'thresh', 0.1);
        
    E_GR = xfm_NUFFT(Mat_size,ones(Mat_size),[],...
            reshape(kspace_GR, [], 1, 1, 2));
    kdata_GR_c = reshape(mean(kdata_GR{S_max_idx,subj},3),[],1,NCoils);
    
    
    for c = 1:NCoils
        disp(['Calculating coil sensitivity map ' num2str(c) ' of ' num2str(NCoils)])
        ims(:,:,:,c) =   fista(E_GR, kdata_GR_c(:,:,c), 1, 0.00000001, ...
            [E_GR.Nd, 1,1], 10, 0.5);
        imshow(abs(ims(:,:,:,c)),[])
        drawnow
    end
    disp(['Adaptive combine estimation of coil sensitivities - GR sampling'])
    sens_GR{subj} = adaptive_estimate_sens('data', permute(ims,[4,1,2,3]),...
            'kernel', 10, 'thresh', 0.1);
        
    E_SILVER = xfm_NUFFT(Mat_size,ones(Mat_size),[],...
            reshape(kspace_SILVER, [], 1, 1, 2));
    kdata_SILVER_c = reshape(mean(kdata_SILVER{S_max_idx,subj},3),[],1,NCoils);
    
    
    for c = 1:NCoils
        disp(['Calculating coil sensitivity map ' num2str(c) ' of ' num2str(NCoils)])
        ims(:,:,:,c) =   fista(E_SILVER, kdata_SILVER_c(:,:,c), 1, 0.00000001, ...
            [E_SILVER.Nd, 1,1], 10, 0.5);
        imshow(abs(ims(:,:,:,c)),[])
        drawnow
    end
    disp(['Adaptive combine estimation of coil sensitivities - SILVER sampling'])
    sens_SILVER{subj} = adaptive_estimate_sens('data', permute(ims,[4,1,2,3]),...
            'kernel', 10, 'thresh', 0.1);
end

save('examples/example7_invivo/example_sensitivities_68_153_306.mat', 'sens_SILVER', 'sens_GR', 'sens_UNIFORM', 'S', '-v7.3');
function [sens] = generate_coil_sensitivity_maps_TURBINE(kdata_filename_base,scaling, tot_spokes,slices)
%%%%%%%%%%%%
% This script takes all 8 .mat slices of TURBINE data (after being read in
% using 'read_in_and_save_TURBINE') and creates slice-wise sensitivity
% maps based on the GR acquisition (1a). scaling can be changed to create
% sensitivity maps with varying spatial resolution (1 means full resolution
% 2mm x 2mm, e.g. 0.5 would mean 4mm x 4mm). It saves both a full
% sensitivity map and one compressed to 8 coils.
%
% Sophie Schauman 2020
%%%%%%%%%%%%%

NFrames = 1;
NSamps = 200;
NCoils = 32;

kspace = gen_radial(roundn(gr2D*180,-3),NSamps, tot_spokes,1,360,1);
kspace = kspace./scaling;
idx = find(abs(kspace(:,1,1))<pi);
kspace =  reshape(kspace(idx,:,:),[],NFrames,2);

E = xfm_NUFFT([100*scaling,100*scaling,1,NFrames],ones(100*scaling),[],kspace, 'wi',1);

    
    
for slice = slices
    
    kdata_filename = [kdata_filename_base num2str(slice)];
    load(kdata_filename);

    kdata = reshape(kdata(idx,:,:,:),[],NFrames,size(kdata,4));
    kdata = kdata.*E.w;% pre-weighting


    for c = 1:size(kdata,3)
        disp(['Calculating coil sensitivity map ' num2str(c) ' of ' num2str(NCoils) ' for slice: ' num2str(slice)])
%         recon(:,:,:,:,c) = E'.*kdata(:,:,c);
        recon(:,:,:,:,c)=reshape(E.iter(reshape(kdata(:,:,c),E.dsize), @pcg, 1E-4, 25),scaling*100,scaling*100,1,1);
        imagesc(abs(recon(:,:,:,:,c)))
        drawnow
    end
    
    
    
    disp(['Adaptive combine estimation of coil sensitivities'])
    sens(:,:,slice,:) = reshape(adaptive_estimate_sens('data', permute(recon,[5,1,2,3,4]),...
            'kernel', round(3*scaling), 'thresh', 0.1),[100*scaling 100*scaling,1,NCoils]);
    
end

psens = compress_coils(sens, 8);

save([kdata_filename(1:end-7) '_sens' num2str(scaling) '.mat'], 'sens', 'psens')


end


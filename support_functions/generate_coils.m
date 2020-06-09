%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulated coil sensitivities
%
%   Sophie Schauman 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear
close all



M = 64; % matrix size
Nc = 32; % Number of coils to generate
Ncc = 8; % Number of coils to compress to
seed = 4; % Random number generator seed


rng(seed)


x = (1:M); y = (1:M)';

for c = 1:Nc
    xcen = randi(M);
    ycen = randi(M);
    sigma = randi([10 M]);

    xcen_p = randi(M);
    ycen_p = randi(M);
    sigma_p = M;

    I_mag = exp(-((x-xcen).^2+(y-ycen).^2)./(2*sigma.^2));
    I_phs = exp(-((x-xcen_p).^2+(y-ycen_p).^2)./(2*sigma_p.^2));
    I_phs = I_phs*4*pi-2*pi;

    psens(:,:,1,c) = I_mag.*exp(1j*I_phs);
end

[u,s,~] = lsvd(reshape(psens,[],Nc),Ncc);
psens = reshape(u*s,M,M,1,Ncc);


save(['examples/example5_gfactor/simu_sensitivity_map' num2str(seed) '.mat'], 'psens')
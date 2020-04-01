% [g, PSF] = gFactor(sens, k)
%
% Computes g-factor noise amplification from sensitivity maps (sens) and
% k-space sampling locations (k)
% Assumes iso-centre is located at floor(N/2)+1
%
% sens - is a [Nx, Ny, Nz, Nc] complex array of coil sensitivities
% k    - is a [Nk,2] or [Nk,3] set of 2D or 3D k-space locations normalised 
%        so that the min/max k-locations in any direction are -1/+1
%
% g    - is a [Nx, Ny, Nz] image of g-factors
% PSF  - is a [Nx*Ny*Nz, Nx*Ny*Nz] point-spread function linear operator

function [g, PSF] = gFactor(sens, k)

   if size(k,2) == 2
        k(:,3)  =   0;
   end
    
    [Nx, Ny, Nz, Nc] = size(sens);
    
    sens    =   reshape(sens,[],Nc);
    idx     =   find(sum(abs(sens).^2,2));
    
    [x,y,z] =   meshgrid(-floor(Nx/2):ceil(Nx/2)-1,...
                         -floor(Ny/2):ceil(Ny/2)-1,...
                         -floor(Nz/2):ceil(Nz/2)-1);
    F       =   exp(-1j*pi*(k(:,1)*x(:)' + k(:,2)*y(:)' + k(:,3)*z(:)'));
    PSF     =   F'*F/size(k,1);
    
    % Theoretically we could do the following, but we don't for memory:
    % g = reshape(diag(inv(S'*blkdiag(PSF,...,PSF)*S)),Nx,Ny,Nz);
    
    tmp     =   zeros(length(idx));
    for c   =   1:Nc
        S   =   spdiags(sens(idx,c),0,length(idx),length(idx));        
        tmp =   tmp + S'*PSF(idx,idx)*S;
    end
    
    g       =   zeros(Nx, Ny, Nz);
    g(idx)  =   sqrt(diag(pinv(tmp)));
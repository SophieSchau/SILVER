function [g, n, varargout] = gFactor_mask(sens, k, mask, ncov)
%
% Computes g-factor noise amplification from sensitivity maps (sens) and
% k-space sampling locations (k). Assuming a reconstruction using the
% pseudo-inverse followed by masking of k-space to suppress signal in 
% unconstrained k-space locations.
% Assumes iso-centre is located at floor(N/2)+1
%
% sens  - is a [Nx, Ny, Nz, Nc] complex array of coil sensitivities
% k     - is a [Nk,2] or [Nk,3] set of 2D or 3D k-space locations normalised 
%         so that the min/max k-locations in any direction are -pi/+pi
% mask  - is mask used in k-space in the final reconstruction step
%
% g     - is a [Nx, Ny, Nz] image of g-factors
% n     - is a [Nx, Ny, Nz] image of predicted noise standard deviation



    if nargin < 4
        ncov = [];
    end
    if nargin < 3
        mask = [];
    end
    if size(k,2) == 2
        k(:,3)  =   0;
    end
    
    [Nx, Ny, Nz, Nc] = size(sens);
    
    sens    =   reshape(sens,[],Nc);
    idx     =   find(sum(abs(sens).^2,2));
    
    [y,x,z] =   meshgrid(-floor(Ny/2):ceil(Ny/2)-1,...
                         -floor(Nx/2):ceil(Nx/2)-1,...
                         -floor(Nz/2):ceil(Nz/2)-1);
    F       =   exp(-1j*(k(:,1)*x(:)' + k(:,2)*y(:)' + k(:,3)*z(:)'));
    PSF_traj     =   F(:,idx)'*F(:,idx);
    clear('F', 'k', 'x', 'y', 'z') % memory management
    
    % Theoretically we could do the following, but we don't for memory:
    % g = reshape(diag(inv(S'*blkdiag(PSF,...,PSF)*S)),Nx,Ny,Nz);
    
    % this is E'E (PSF including coils):
    PSF_traj_coils = conj(sens(idx,:)*sens(idx,:)').*PSF_traj; 
    
    if ~isempty(ncov)
        tmp= 0; % This is E'nn'E <- noise covariance included
        for i = 1:Nc
            tmp = tmp + conj(sens(idx,i)).*(PSF_traj.*(ncov(i,:)*sens(idx,:).'));
        end
    else
        tmp = PSF_traj_coils;
    end
    clear('PSF_traj') % memory management
    
    if ~isempty(mask)
        M = zeros(Nx*Ny*Nz); % convolution operator (because the mask is point-wise multiplied in k-space)
        n = 0;
        % get M, by generating it column by column
        for i = 1:Ny
            for j = 1:Nx
                n = n+1;
                tmp_im = zeros(Nx,Ny);
                tmp_im(j,i) = 1;
                M(:,n) = reshape(ifft2(ifftshift(fftshift(fft2(tmp_im)).*mask)),[],1);
            end
        end
    else
        M = 1;
    end
    
    n       =   zeros(Nx, Ny, Nz);
    g       =   zeros(Nx, Ny, Nz);
    
    
    tmp2 = zeros(prod([Nx,Ny,Nz]));
    tmp2(idx,idx) = pinv(PSF_traj_coils)*tmp*pinv(PSF_traj_coils); % insert noise amplification calculation without masking into fully sized matrix
    
    n  =   reshape(sqrt(diag(M*tmp2*M')),Nx,Ny,Nz); % do left and right multioplication with mask convolution operator and take sqrt of the diagonal for noise variance predictions
    
    %%%%%%%%% not sure this is right...
    if isempty(ncov)
        g(idx)  =   n(idx).*sqrt(diag(PSF_traj_coils));
    else
        g(idx)  =   sqrt(Nx*Ny)*n(idx)./sqrt(diag(sens(idx,:)*ncov*sens(idx,:)')./sum(abs(sens(idx,:)).^2,2));
    end
    %%%%%%%%%%%
    

end
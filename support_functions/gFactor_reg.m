% [g, n] = gFactor_reg(sens, k, lambda)
%
% Computes g-factor noise amplification from sensitivity maps (sens) and
% k-space sampling locations (k)
% Assumes iso-centre is located at floor(N/2)+1
%
% sens  - is a [Nx, Ny, Nz, Nc] complex array of coil sensitivities
% k     - is a [Nk,2] or [Nk,3] set of 2D or 3D k-space locations normalised 
%         so that the min/max k-locations in any direction are -pi/+pi
% lambda- is the scalar reguarlisation parameter
%
% g     - is a [Nx, Ny, Nz] image of g-factors
% n     - is a [Nx, Ny, Nz] image of predicted noise standard deviation

function [g, n, varargout] = gFactor_reg(sens, k, lambda, ncov)

    if nargin < 4
        ncov = [];
    end
    if nargin < 3
        lambda=0;
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
    PSF     =   F(:,idx)'*F(:,idx);
    clear('F', 'k', 'x', 'y', 'z') % memory management
    
    % Theoretically we could do the following, but we don't for memory:
    % g = reshape(diag(inv(S'*blkdiag(PSF,...,PSF)*S)),Nx,Ny,Nz);
    
    tmp = conj(sens(idx,:)*sens(idx,:)').*PSF;
    if ~isempty(ncov)
        tmp2= 0;
        for i = 1:Nc
            %tmp3 = F(:,idx).*(ncov(i,:)*sens(idx,:).');
            tmp2 = tmp2 + conj(sens(idx,i)).*(PSF.*(ncov(i,:)*sens(idx,:).'));
        end
    else
        tmp2 = tmp;
    end
    clear('PSF')
    n       =   zeros(Nx, Ny, Nz);
    g       =   zeros(Nx, Ny, Nz);
    n(idx)  =   sqrt(diag((tmp + lambda*eye(length(idx)))\tmp2/(tmp + lambda*eye(length(idx)))));
    if isempty(ncov)
        g(idx)  =   n(idx).*sqrt(diag(tmp2));
    else
        g(idx)  =   sqrt(Nx*Ny)*n(idx)./sqrt(diag(sens(idx,:)*ncov*sens(idx,:)')./sum(abs(sens(idx,:)).^2,2));
    end
    
    

end
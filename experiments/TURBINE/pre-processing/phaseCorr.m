function [data varargout] = phaseCorr(pc_data, data, varargin)

%   Mark Chiew
%   September 2014
%
%   Assumes squeezed data format that is output from mapVBVD
%
%   [x,     (1)
%    coils, (2)
%    y,     (3)
%    par,   (4) 
%    sli,   (5)
%    ave,   (6)
%    phs,   (7)
%    eco,   (8)
%    rep,   (9)
%    set,   (10)
%    seg    (11)]


p   =   inputParser;

p.addParameter('x',     0.5,    @isscalar);
p.addParameter('mode', '3D');

p.parse(varargin{:});

x       =   p.Results.x;
mode    =   p.Results.mode;

nc      =   size(pc_data,2);
nx      =   size(pc_data,1);
ny      =   size(data,3);
nz      =   size(data,4);
nMeas   =   size(pc_data,9);
p       =   zeros(nc,2);
w       =   zeros(nc,1);
pcorr   =   [];
ns      =   size(data,10);


switch mode
case 'CAIPI'
    pc_data =   sum(pc_data, 4);
    for m = 1:nMeas
        for z = 1:ns
            Sp1 =   ifftdim(squeeze(pc_data(:,:,1,1,1,1,1,1,m,z,2)),1);
            Sp2 =   ifftdim(squeeze(pc_data(:,:,1,1,1,2,1,1,m,z,2)),1);
            Sp  =   0.5*(Sp1+Sp2);
            Sn  =   ifftdim(squeeze(pc_data(:,:,1,1,1,1,1,1,m,z,1)),1);

            S       =   Sn.*conj(Sp);
            S       =   mean(S,2);
            
            x0          =   find(abs(sqrt(S(:)))./max(sqrt(abs(S(:)))) > 0.1);
            x0          =   (x0(1):x0(end))';
            phi         =   unwrap(angle(S(x0)));
            [p,stat]    =   polyfit(x0, phi, 1);

            pcorr(z, m)   =   p(1);

            phi =   polyval(p, 1:nx).';
            if nargout > 2
                ref(:,:,1,z,1,1,1,1,m,1,2)  =   fftdim(Sn.*exp(-x*j*repmat(phi,1,nc)),1);
                ref(:,:,1,z,1,1,1,1,m,1,1)  =   fftdim(Sp1.*exp((1-x)*j*repmat(phi,1,nc)),1);
                ref(:,:,1,z,1,1,1,1,m,1,3)  =   fftdim(Sp2.*exp((1-x)*j*repmat(phi,1,nc)),1);
            end
            if ~isempty(data)
                data(:,:,:,:,1,1,1,1,m,z,1) =   bsxfun(@times, ifftdim(data(:,:,:,:,1,1,1,1,m,z,1),1), exp(-x*j*phi));
                data(:,:,:,:,1,1,1,1,m,z,2) =   bsxfun(@times, ifftdim(data(:,:,:,:,1,1,1,1,m,z,2),1), exp((1-x)*j*phi));
            end
        end
    end
case 'TURBINE'

    for m = 1:nMeas
    for z = 1:nz
        Sp1 =   ifftdim(squeeze(pc_data(:,:,1,z,1,1,1,1,m,1,2)),1);
        Sp2 =   ifftdim(squeeze(pc_data(:,:,1,z,1,2,1,1,m,1,2)),1);
        Sp  =   0.5*(Sp1+Sp2);
        Sn  =   ifftdim(squeeze(pc_data(:,:,1,z,1,1,1,1,m,1,1)),1);

        S   =   Sn.*conj(Sp);
        S   =   mean(S,2);
        
        x0          =   find(abs(sqrt(S(:)))./max(abs(sqrt(S(:)))) > 0.1);
        x0          =   (x0(1):x0(end))';
        phi         =   unwrap(angle(S(x0)));
        [p,stat]    =   polyfit(x0, phi, 1);

        pcorr(z,m)  =   p(1);
        
        phi =   polyval(p, 1:nx).';

        if nargout > 2
            ref(:,:,1,z,1,1,1,1,m,1,2)  =   fftdim(bsxfun(@times, Sn, exp(-x*j*phi)),1);
            ref(:,:,1,z,1,1,1,1,m,1,1)  =   fftdim(bsxfun(@times, Sp1, exp((1-x)*j*phi)),1);
            ref(:,:,1,z,1,1,1,1,m,1,3)  =   fftdim(bsxfun(@times, Sp2, exp((1-x)*j*phi)),1);
        end
        if ~isempty(data)
            data(:,:,:,z,1,1,1,1,m,1,1) =   fftdim(bsxfun(@times, ifftdim(data(:,:,:,z,1,1,1,1,m,1,1),1), exp(-x*j*phi)),1);
            data(:,:,:,z,1,1,1,1,m,1,2) =   fftdim(bsxfun(@times, ifftdim(data(:,:,:,z,1,1,1,1,m,1,2),1), exp((1-x)*j*phi)),1);
        end
    end
    end

case '3D'
    pc_data =   sum(pc_data,4);
    pc_data =   sum(pc_data, 10);
    if size(pc_data,9) == 1

        Sp1 =   ifftdim(squeeze(pc_data(:,:,1,1,1,1,1,1,1,1,2)),1);
        Sp2 =   ifftdim(squeeze(pc_data(:,:,1,1,1,2,1,1,1,1,2)),1);
        Sp  =   0.5*(Sp1+Sp2);
        Sn  =   ifftdim(squeeze(pc_data(:,:,1,1,1,1,1,1,1,1,1)),1);

        S   =   Sn.*conj(Sp);
        S   =   mean(S,2);
            
        
        x0          =   find(abs(S(:))./max(abs(S(:))) > 0.1);
        phi         =   unwrap(angle(S(x0)));
        [p,stat]    =   polyfit(x0, phi, 1);

        pcorr   =   p(1);

        phi =   polyval(p, 1:nx).';

        if nargout > 2
            ref(:,:,1,:,1,1,1,1,:,1,2)  =   fftdim(bsxfun(@times, Sn, exp(-x*j*phi)),1);
            ref(:,:,1,:,1,1,1,1,:,1,1)  =   fftdim(bsxfun(@times, Sp1, exp((1-x)*j)),1);
            ref(:,:,1,:,1,1,1,1,:,1,3)  =   fftdim(bsxfun(@times, Sp2, exp((1-x)*j)),1);
        end
        if ~isempty(data)
            data(:,:,:,:,1,1,1,1,:,1,1) =   fftdim(bsxfun(@times, ifftdim(data(:,:,:,:,1,1,1,1,:,1,1),1), exp(-x*j*phi)),1);
            data(:,:,:,:,1,1,1,1,:,1,2) =   fftdim(bsxfun(@times, ifftdim(data(:,:,:,:,1,1,1,1,:,1,2),1), exp((1-x)*j*phi)),1);
        end

    else

        for m = 1:nMeas
            Sp1 =   ifftdim(squeeze(pc_data(:,:,1,1,1,1,1,1,m,1,2)),1);
            Sp2 =   ifftdim(squeeze(pc_data(:,:,1,1,1,2,1,1,m,1,2)),1);
            Sp  =   0.5*(Sp1+Sp2);
            Sn  =   ifftdim(squeeze(pc_data(:,:,1,1,1,1,1,1,m,1,1)),1);

            S   =   Sn.*conj(Sp);
            %[S,~,~]=    lsvd(S,1);
            S   =   mean(S,2);
            
            x0          =   find(abs(S(:))./max(abs(S(:))) > 0.1);
            phi         =   unwrap(angle(S(x0)));
            [p,stat]    =   polyfit(x0, phi, 1);

            pcorr   =   p(1);

            phi =   polyval(p, 1:nx).';

            if nargout > 2
                ref(:,:,1,:,1,1,1,1,m,1,2)  =   fftdim(bsxfun(@times, Sn, exp(-x*j*phi)),1);
                ref(:,:,1,:,1,1,1,1,m,1,1)  =   fftdim(bsxfun(@times, Sp1, exp((1-x)*j)),1);
                ref(:,:,1,:,1,1,1,1,m,1,3)  =   fftdim(bsxfun(@times, Sp2, exp((1-x)*j)),1);
            end
            if ~isempty(data)
                data(:,:,:,:,1,1,1,1,m,1,1) =   fftdim(bsxfun(@times, ifftdim(data(:,:,:,:,1,1,1,1,m,1,1),1), exp(-x*j*phi)),1);
                data(:,:,:,:,1,1,1,1,m,1,2) =   fftdim(bsxfun(@times, ifftdim(data(:,:,:,:,1,1,1,1,m,1,2),1), exp((1-x)*j*phi)),1);
            end
        end
    end
case '2D'
    nz = size(pc_data, 5);
    for m = 1:nMeas
        for z = 1:nz
            Sp1 =   ifftdim(squeeze(pc_data(:,:,1,1,z,1,1,1,m,1,2)),1);
            Sp2 =   ifftdim(squeeze(pc_data(:,:,1,1,z,2,1,1,m,1,2)),1);
            Sp  =   0.5*(Sp1+Sp2);
            Sn  =   ifftdim(squeeze(pc_data(:,:,1,1,z,1,1,1,m,1,1)),1);

            S   =   Sn.*conj(Sp);
            S   =   mean(S,2);
            
            x0          =   find(abs(sqrt(S(:)))./max(sqrt(abs(S(:)))) > 0.1);
            x0          =   (x0(1):x0(end))';
            phi         =   unwrap(angle(S(x0)));
            [p,stat]    =   polyfit(x0, phi, 1);

            pcorr(z, m)   =   p(1);
            
            phi =   polyval(p, (1:nx)).';
            corr=   repmat(phi,[1,nc, ny]);
            if nargout > 2
                ref(:,:,1,1,z,1,1,1,m,1,2)  =   fftdim(Sn.*exp(-x*j*repmat(phi,1,nc)),1);
                ref(:,:,1,1,z,1,1,1,m,1,1)  =   fftdim(Sp1.*exp((1-x)*j*repmat(phi,1,nc)),1);
                ref(:,:,1,1,z,1,1,1,m,1,3)  =   fftdim(Sp2.*exp((1-x)*j*repmat(phi,1,nc)),1);
            end
            if ~isempty(data)
                data2   =   fftdim(ifftdim(data(:,:,:,1,z,1,1,1,m,1,1),1).*exp(-x*j*corr),1);
                data(:,:,:,1,z,1,1,1,m,1,1) =   data2;
                data2   =   fftdim(ifftdim(data(:,:,:,1,z,1,1,1,m,1,2),1).*exp((1-x)*j*corr),1);
                data(:,:,:,1,z,1,1,1,m,1,2) =   data2;

            end
        end
    end
end

data    =   sum(data,11);
if nargout > 1
    varargout{1}    =   pcorr;
end

if nargout > 2
    varargout{2}    =   squeeze(ref);
end

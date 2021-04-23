function [data varargout] = physioCorr(data,ref,order)

%   Mark Chiew
%   Feb 2015
%

if nargin < 2
    ref =   [];
end
if nargin < 3
    order = 0;
end

nc      =   size(data,2);
nx      =   size(data,1);
ny      =   size(data,3);
nt      =   size(data,4);

% Get y-indices
y   =   (1:ny)'-floor(ny/2)+1;      

% iFFT along z-direction
%data=   ifftdim(data,3);

% iFFT temporal mean as reference, along x-direction
% and get centre of k-space column

if isempty(ref)
    ref  =  mean(data,4); 
end
ref  =   sum(ifftdim(ref,1).*repmat(exp(+1j*2*pi*0.5*((1:nx)-(nx/2+1))'/nx),[1,nc,ny]),1);
%{
if isempty(ref)
    ref  =  1;
end
%}

% Loop over all blades
for i = 1:nt
    % iFFT first blade as reference, along x-direction
    % and get centre of k-space column
    d   =   ifftdim(data(:,:,:,i),1);
    d   =   sum(d.*repmat(exp(+1j*2*pi*0.5*((1:nx)-(nx/2+1))'/nx),[1,nc,ny]),1);

    %{
    % mean phase difference across all coils
    d   =   squeeze(mean(d.*conj(ref),2));

    % find voxels with signal
    yy  =   find(sqrt(abs(d))>0.25*max(sqrt(abs(d))));
    yy  =   (yy(1):yy(end))';

    % fit linear phase model
    phi =   unwrap(angle(d(yy)));
    p   =   polyfit(y(yy)-(floor(ny/2)+1), phi, max(order,0));

    % save fit coefficients
    pcorr(:,i) =   double(p);
    %}
    d   =   squeeze(mean(mean(d.*conj(ref),2),3));
    pcorr(:,i)  =   double(angle(d));
   
end
%{
[b, a]      =   butter(10, 0.05, 'low');
pcorr(1,:)  =   filtfilt(b, a, unwrap(pcorr(1,:)));
pcorr(2,:)  =   filtfilt(b, a, unwrap(pcorr(2,:)));
%}

% Apply phase correction
for i = 1:nt
    corr    =   reshape(polyval(pcorr(:,i), y),1,1,[]);
    corr    =   repmat(exp(-1j*corr),nx,nc,1);
    data(:,:,:,i)=data(:,:,:,i).*corr;
end

if nargout > 1
    varargout{1}    =   pcorr;
end

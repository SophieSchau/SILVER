function Phis = gr3D(idx)

if nargin <1
    idx = 1:2;
end

% Define the Fibonacci matrix as per Chan et al MRM 2009
M = [0 1 0; ...
     0 0 1; ...
     1 0 1];
 
% Derive the real eigenvector
[V,~] = eig(M);

% Normalisation is different in Chan paper, so force last component to be 1
Phis = V(1:2,1)./V(3,1);

if idx==1
    Phis=Phis(1);
elseif idx==2
    Phis=Phis(2);
elseif idx == 1:2
    Phis = Phis;
else
    error('Only 1, and 2 are valid inputs to this funcion')
end
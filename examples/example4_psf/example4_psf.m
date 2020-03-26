%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   What do the point-spread-functions of trajectories optimized with
%   SILVER look like compared with golden ratio method or uniform radial
%   sampling?
%                                              
%   Sophie Schauman 2020                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% 1. Choose  set of window sizes, S, to consider
S = [16,32,48];

%% 2. Do SILVER optimization for that range
for s = S
   savename = ['examples/precalculated/silver_' strrep(num2str(s{:}),' ', '_') '.mat'];
   if ~exist(savename, 'file')
        SILVER_2D(s{:},'electrostatic_potential',savename) ;
   end
   disp([num2str(s{:}) ' done'])
end

%% 3. Generate PSF's for SILVER, GR, and uniform sampling

for s = S
    
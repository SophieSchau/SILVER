function [energy, coordinates] = low_energy_configuration(N)
%LOW_ENERGY_CONFIGURATION Numerical optimisation of N spokes using
%electrostatic potential by treatng thw spoke tips as unit charges on the
%unit sphere
%
%   INPUTS:   N           - number of spokes
%   OUTPUTS:  energy      - the lowest potential energy found
%             coordinates - the configuration that gave the lowest energy
%
%   DEPENDENCIES: 
%             electrostatic_potential()
%
%
% Sophie Schauman, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, check if this has been optimised already
if exist([num2str(N) 'spokes_config.mat'], 'file')
    load([num2str(N) 'spokes_config.mat'], 'energy', 'coordinates')
else
    
starting_coordinates = rand(N,2);
options = optimset('Display','iter', 'MaxFunEvals', 100000);
[coordinates, energy] = fmincon(@electrostatic_potential,starting_coordinates, [], [], [], [], zeros(size(starting_coordinates)), ones(size(starting_coordinates)), [], options); 
 
try
    save(['efficiency_estimation_3D/low_energy_configurations/' num2str(N) 'spokes_config.mat'], 'energy', 'coordinates')
catch
    save([num2str(N) 'spokes_config.mat'], 'energy', 'coordinates')
end

end
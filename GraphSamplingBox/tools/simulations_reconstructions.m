function simulations_reconstructions(graph, noise, reg_term)
% Script to launch the reconstruction simulations
%
% Copyright (c) 2016 G. Puy, N. Tremblay
%
% This file is part of the GraphSamplingBox
%
% The GraphSamplingBox is free software: you can redistribute it and/or 
% modify it under the terms of the GNU Affero General Public License as 
% published by the Free Software Foundation, either version 3 of the 
% License, or (at your option) any later version.
%
% The GraphSamplingBox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% If you use this toolbox please kindly cite
%  [1] G. Puy, N. Tremblay, R. Gribonval and P. Vandergheynst. "Random 
%   sampling of bandlimited signals on graphs", ACHA, 2016.

% --- Parameters for simulations
param_sim.band_limit = 10; % Band limit
param_sim.nb_meas = 200; % Nb of measurements
param_sim.nb_sim = 10; % Nb of simul. per set of parameters
param_sim.reg = logspace(-3, 2, 50); % Default set of regularisation parameters

% --- Simulations with or without noise 
if strcmp(noise, 'nonoise')
    param_sim.noise_list = 0;
else
    param_sim.noise_list = logspace(-3, -2+log10(5), 10);
end

% --- Set regularisation term
param_sim.reg_term = reg_term;

% --- Construct graph
switch graph
    case 'minnesota'
        G = gsp_minnesota;
        % Specific set of regularisation parameters for minnesota
        param_sim.reg = logspace(-1, 10, 50);
    case 'bunny'
        G = gsp_bunny;
    case 'community'
        % Graph size
        N = 1000;
        % Nb. of clusters
        param_graph.Nc = param_sim.band_limit;
        % Size of clusters
        param_graph.com_sizes = round(N/param_graph.Nc/10);
        param_graph.com_sizes = [param_graph.com_sizes, ...
            repmat(floor((N-param_graph.com_sizes)/(param_graph.Nc-1)), 1, ...
            param_graph.Nc-2)];
        param_graph.com_sizes = [param_graph.com_sizes, ...
            N-sum(param_graph.com_sizes)];
        param_graph.com_lims = [0, cumsum(param_graph.com_sizes)];
        % Probability of connections intercluster
        param_graph.world_density = 1/N;
        % Construct graph
        G = gsp_community(N, param_graph);    
end


% --- Launch simulations
[err, errUk, errBarUk, timeRecons] = ...
    launch_simulations_of_reconstructions(G, param_sim);

% --- Save results
save(['results/rec_', graph, '_', reg_term, '_', noise, '.mat'], ...
    'G', 'err', 'errUk', 'errBarUk', 'timeRecons', 'param_sim');

end
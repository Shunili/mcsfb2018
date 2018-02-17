function simulations_eigenvalues(graph, band_limit, sampling, factor)
% Script to prepare the eigenvalue simulations to test if the RIP is
% satisfied.
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
param_sim.band_limit = band_limit; % Band limit
param_sim.nb_sim = 500; % Nb of simul. per set of parameters

% --- Construct graph
switch graph
    case 'minnesota'
        G = gsp_minnesota;
    case 'bunny'
        G = gsp_bunny;
    case 'comet'
        G = gsp_comet(500, 495);
    case 'tree'
        G = gsp_tree(9,2);
    case 'path'
        G = gsp_path(1000);
    case 'community'
        % Graph size
        N = 1000;
        % Nb. of clusters
        param_graph.Nc = param_sim.band_limit;
        % Size of clusters
        param_graph.com_sizes = round(N/param_graph.Nc/factor);
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

% --- Number of measurements
if band_limit < 50
    param_sim.nb_meas = band_limit:2:600;
elseif band_limit >= 50
    param_sim.nb_meas = band_limit:50:2000;
end

% --- Launch simulations
param_sim.sampling = sampling;
[max_eig, min_eig] = launch_simulations_eigenvalues(G, param_sim);

% --- Save results
G = gsp_compute_fourier_basis(G);
if strcmp(graph, 'community')
    save(['results/eig_', graph, '_', int2str(band_limit), '_', ...
        sampling, '_', int2str(factor), '.mat'], ...
        'param_sim', 'max_eig', 'min_eig', 'G');
else
    save(['results/eig_', graph, '_', int2str(band_limit), '_', ...
        sampling, '.mat'], 'param_sim', 'max_eig', 'min_eig', 'G');
end

end
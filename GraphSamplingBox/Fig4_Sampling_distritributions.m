% This script allows one to reproduce Fig. 4 of [1]
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

clc;
clear;
setup;

%% Options
% Choose the graph: 'minnesota', 'bunny' or 'community'
graph = 'minnesota';
% Choose the band limit
band_limit = 10;
% Choose 'optimal' for optimal sampling or 'estimated' for estimation with
% Algorithm 1 of [1]
sampling = 'estimated'; 

%% Construct graph
switch graph
    case 'minnesota'
        G = gsp_minnesota;
        cp = [-93.49 45 12];
    case 'bunny'
        G = gsp_bunny;
        cp = [0.13 0 20];
    case 'community'
        % Graph size
        N = 1000;
        % Parameter controlling the size of the smallest community
        factor = 8;
        % Nb. of clusters
        param_graph.Nc = band_limit;
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
        cp = [-2.9328, -2.5947, 17.3205];
end

%% Sampling
switch sampling
    case 'estimated'
        param_est.nb_estimation = 1;
        [~, cumCoh] = gsp_estimate_lk(G, band_limit, ...
            param_est);
        weight = cumCoh/sum(cumCoh(:));
    case 'optimal'
        G = gsp_compute_fourier_basis(G);
        Uk = G.U(:, 1:band_limit);
        weight = sum(Uk.^2, 2)/sum(Uk(:).^2);
end

%% Show result
figure(1); clf;
param.show_edges = 0;
gsp_plot_signal(G, weight, param);
colormap(flipud(hot));
set(gca,'CameraPosition', cp)
set(gcf, 'color', 0.8*ones(1, 3),'inverthardcopy', 'off');
axis square;
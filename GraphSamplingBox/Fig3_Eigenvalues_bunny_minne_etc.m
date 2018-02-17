% This script allows one to reproduce Fig. 3 of [1]. 
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

%%
% Choose the graph: 'minnesota', 'bunny', 'path'
graph = 'minnesota';
list_sampling = {'uniform', 'optimal', 'estimated'};
line_color = {'k', 'r', 'b'};

%%
list_option = {10, 100}; % Band-limits

%% Run the simulations (comment this part if results are already available)
for ind_s = 1:numel(list_sampling)
    for ind_f = 1:numel(list_option)
        simulations_eigenvalues(graph, list_option{ind_f}, ...
            list_sampling{ind_s});
    end
end

%% Load results
%
meas = cell(numel(list_sampling), numel(list_option));
prob = cell(numel(list_sampling), numel(list_option));
%
for ind1 = 1:numel(list_sampling)
    for ind2 = 1:numel(list_option)
        %
        filename = ['results/eig_', graph, '_', ...
            int2str(list_option{ind2}), '_', ...
            list_sampling{ind1}, '.mat'];
        %
        if exist(filename, 'file')
            load(filename);
        else
            continue;
        end
        %
        prob{ind1, ind2} = sum(min_eig>0.005)/param_sim.nb_sim;
        meas{ind1, ind2} = param_sim.nb_meas;
    end
end

%% Show results
line_style = {'-+', '-+', '-+'};
lim_x = [200 2000];
%
for ind2 = 1:numel(list_option)
    
    figure(ind2); clf; hold on;
    set(gca, 'ytick', 0:.2:1);
    set(gca, 'fontsize', 16);
    xlabel('m', 'fontsize', 16); ylabel('1 - \epsilon', 'fontsize', 16);
    xlim([param_sim.band_limit, lim_x(ind2)]); 
    ylim([0 1]); axis square;
    
    for ind1 = 1:numel(list_sampling)
        
        plot(meas{ind1, ind2}, prob{ind1, ind2}, line_style{ind2}, ...
            'Color', line_color{ind1}, 'linewidth', 1.5);
        
    end
end
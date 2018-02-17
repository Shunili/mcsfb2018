% This script allows one to reproduce Fig. 2 of [1]. 
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
graph = 'community';
list_sampling = {'uniform', 'optimal', 'estimated'};

%%
list_option = {1, 2, 4, 6, 8};

%% Run the simulations (comment this part if results are already available)
for ind_s = 1:numel(list_sampling)
    for ind_f = 1:numel(list_option)
        simulations_eigenvalues(graph, 10, list_sampling{ind_s}, ...
            list_option{ind_f});
    end
end

%% Load results
%
cumCoh = cell(numel(list_sampling), numel(list_option));
meas = cell(numel(list_sampling), numel(list_option));
prob = cell(numel(list_sampling), numel(list_option));
%
for ind1 = 1:numel(list_sampling)
    for ind2 = 1:numel(list_option)
        %
        filename = ['results/eig_', graph, '_10_', ...
            list_sampling{ind1}, '_', ...
            int2str(list_option{ind2}), '.mat'];
        %
        if exist(filename, 'file')
            load(filename);
        else
            continue;
        end
        %
        prob{ind1, ind2} = sum(min_eig>0.005)/param_sim.nb_sim;
        meas{ind1, ind2} = param_sim.nb_meas;
        %
        switch list_sampling{ind1}
            case 'uniform'
                cumCoh{ind1, ind2} = G.N*max(sum(G.U(:, 1:10).^2, 2));
            case 'optimal'
                cumCoh{ind1, ind2} = 10;
            case 'estimated'
                cumCoh{ind1, ind2} = 10;
        end
    end
end

%% Show results
line_style = {'-+', '-+', '-+', '-+', '-+'};
line_color = {'k', 'r', 'b', 'g', [1 .5 0]};
%
for ind1 = 1:numel(list_sampling)
    
    figure(ind1); clf; hold on;
    set(gca, 'ytick', 0:.2:1);
    set(gca, 'fontsize', 14);
    xlabel('m', 'fontsize', 14); ylabel('1 - \epsilon', 'fontsize', 14);
    xlim([10 400]); ylim([0 1]); axis square;
    
    for ind2 = 1:numel(list_option)
        plot(meas{ind1, ind2}, prob{ind1, ind2}, line_style{ind2}, ...
            'Color', line_color{ind2}, 'linewidth', 1.5);
        plot(3*cumCoh{ind1, ind2}*[1, 1], [0 1], '--', ...
            'Color', line_color{ind2}, 'linewidth', 1.5);
    end
end
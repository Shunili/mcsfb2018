% This script allows one to reproduce Fig. 5 of [1]
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
% Choose the graph: 'minnesota', 'bunny', 'community'
graph = 'bunny';
noise = 'nonoise';
list_reg = {'L', 'L2', 'L4'};

%% Run the simulations (comment this part if results are already available)
for ind_r = 1:numel(list_reg)
        simulations_reconstructions(graph, noise, list_reg{ind_r})
end

%% Load results
%
err_all = [];
errUk_all = [];
errBarUk_all = [];
timeRecons_all = [];
%
for ind = 1:numel(list_reg)
    %
    filename = ['results/rec_', graph, '_', list_reg{ind}, '_', ...
        noise, '.mat'];
    if exist(filename, 'file')
        load(filename)
    else
        continue;
    end
    err_all = [err_all; mean(squeeze(err))];
    errUk_all = [errUk_all; mean(squeeze(errUk))];
    errBarUk_all = [errBarUk_all; mean(squeeze(errBarUk))];
    timeRecons_all = [timeRecons_all; mean(squeeze(timeRecons))];
end

%%
colors = {'k', 'b', 'r'};
%
figure(1); clf;
for ind = 1:3
    loglog(param_sim.reg, err_all(ind, :), '+-', ...
        'color', colors{ind}, 'linewidth', 1.5); hold on;
end
xlabel('\gamma', 'fontsize', 14);
ylabel('|| x - x* ||_2', 'fontsize', 14);
%
figure(2); clf;
for ind = 1:3
    loglog(param_sim.reg, errUk_all(ind, :), '+-', ...
        'color', colors{ind}, 'linewidth', 1.5); hold on;
end
xlabel('\gamma', 'fontsize', 14);
ylabel('|| \alpha* - x ||_2', 'fontsize', 14);
%
figure(3); clf;
for ind = 1:3
    loglog(param_sim.reg, errBarUk_all(ind, :), '+-', ...
        'color', colors{ind}, 'linewidth', 1.5); hold on;
end
xlabel('\gamma', 'fontsize', 14);
ylabel('|| \beta* ||_2', 'fontsize', 14);
%
figure(4); clf;
for ind = 1:3
    loglog(param_sim.reg, timeRecons_all(ind, :), '+-', ...
        'color', colors{ind}, 'linewidth', 1.5); hold on;
end
for ind = 1:4
    figure(ind);
    set(gca, 'fontsize', 14);
    xlim([param_sim.reg(1), param_sim.reg(end)]);
    ylim([3.5e-5, 1]);
    axis square;
end
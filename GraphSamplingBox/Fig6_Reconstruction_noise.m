% This script allows one to reproduce Fig. 6 of [1]
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
% Choose the graph: 'community' or 'bunny'
graph = 'community';
noise = 'noise';
reg = 'L4';

%% Run the simulations (comment this part if results are already available)
simulations_reconstructions(graph, noise, reg)
% NOTE: To obtain the curve corresponding to the minnesota graph, one
% should manually edit the function simulations_reconstructions in the
% folder tools and change the line 29 from:
%   param_sim.reg = logspace(-3, 2, 50);
% to
%   param_sim.reg = logspace(-1, 10, 50);

%% Load results
%
filename = ['results/rec_', graph, '_', reg, '_', ...
    noise, '.mat'];
if exist(filename, 'file')
    load(filename)
else
    return;
end
%
err_all = squeeze(mean(err));
errUk_all = squeeze(mean(errUk));
errBarUk_all = squeeze(mean(errBarUk));
timeRecons_all = squeeze(mean(timeRecons));
[min_err, I] = min(err_all, [], 1);

%%
colors = {'b', 'r', 'k', 'g', 'c'};
%
figure(1); clf;
for ind = 2:2:10
    loglog(param_sim.reg, err_all(:, ind), '+-', ...
        'color', colors{ind/2}, 'linewidth', 1.5); hold on;
end
loglog(param_sim.reg(I(2:2:10)), min_err(2:2:10), '-o', ...
    'color', [1 .5 0], 'linewidth', 1.5);
xlabel('\gamma', 'fontsize', 14);
ylabel('|| x - x* ||_2', 'fontsize', 14);
%
figure(2); clf;
for ind = 2:2:10
    loglog(param_sim.reg, errUk_all(:, ind), '+-', ...
        'color', colors{ind/2}, 'linewidth', 1.5); hold on;
end
%
figure(3); clf;
for ind = 2:2:10
    loglog(param_sim.reg, errBarUk_all(:, ind), '+-', ...
        'color', colors{ind/2}, 'linewidth', 1.5); hold on;
end
for ind = 1:3
    figure(ind);
    set(gca, 'fontsize', 14);
    xlim([param_sim.reg(1), param_sim.reg(end)]);
    ylim([1e-2, 3e1]);
    axis square;
end
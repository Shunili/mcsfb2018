function [max_eig, min_eig] = launch_simulations_eigenvalues(G, param)
% Script to launch the eigenvalue simulations to test if the RIP is
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


% Parameters
if nargin < 2
    param = struct;
end
if ~isfield(param, 'nb_sim')
    param.nb_sim = 1000;
end
if ~isfield(param, 'band_limit')
    param.band_limit = 10;
end

% Maximum eigenvalue
G = gsp_estimate_lmax(G);

% Fourier basis
G = gsp_compute_fourier_basis(G);
Uk = G.U(:, 1:param.band_limit);

% Init.
min_eig = zeros(param.nb_sim, numel(param.nb_meas));
max_eig = zeros(param.nb_sim, numel(param.nb_meas));

% Simulations
for ind_sim = 1:param.nb_sim

    % Sampling procedure
    switch param.sampling
        case 'uniform'
            weight = ones(G.N, 1)/G.N;
        case 'optimal'
            weight = sum(Uk.^2, 2)/sum(Uk(:).^2);
        case 'estimated'
            param_est.nb_estimation = 1;
            [~, cumCoh] = gsp_estimate_lk(G, param.band_limit, ...
                        param_est);
            weight = cumCoh/sum(cumCoh(:));
    end
    P = sparse(1:G.N, 1:G.N, 1./sqrt(weight), G.N, G.N);

    % For each noise level
    for ind_meas = 1:numel(param.nb_meas)

        % Sampling matrix
        ind_obs = datasample(1:G.N, param.nb_meas(ind_meas), ...
            'Replace', true, 'Weights', weight);
        M = sparse(1:param.nb_meas(ind_meas), ind_obs, 1, ...
            param.nb_meas(ind_meas), G.N);
        
        % Reconstruct signal for all regularisation parameters
        A = M*P*Uk;
        A = (A'*A)/param.nb_meas(ind_meas);
        e = eig(A);
        max_eig(ind_sim, ind_meas) = max(e);
        min_eig(ind_sim, ind_meas) = min(e);
        
    end
end


end

function [err, errUk, errBarUk, timeRecons] = ...
    launch_simulations_of_reconstructions(G, param)
% Script to prepare the reconstruction simulations
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
    param.nb_sim = 10;
end
if ~isfield(param, 'reg')
    param.reg = logspace(-3, 2, 50);
end
if ~isfield(param, 'band_limit')
    param.band_limit = 10;
end
if ~isfield(param, 'nb_meas')
    param.nb_meas = 5*param.band_limit;
end
if ~isfield(param, 'reg_term')
    param.reg_term = 'L';
end
if ~isfield(param, 'noise_list')
    param.noise_list = 0;
end

% Maximum eigenvalue
G = gsp_estimate_lmax(G);

% Fourier basis
G = gsp_compute_fourier_basis(G);
Uk = G.U(:, 1:param.band_limit);

% Init.
err = zeros(param.nb_sim, numel(param.reg), numel(param.noise_list));
errUk = err; errBarUk = err; timeRecons = err;


% Simulations
for ind_sim = 1:param.nb_sim
    
    % For each noise level
    for ind_noise = 1:numel(param.noise_list)
        
        % Estimate lambda_k
        param_est.nb_estimation = 1;
        [lk, cumCoh] = gsp_estimate_lk(G, param.band_limit, param_est);
        
        % Sim. with optimal sampling (WITH replacement)
        weight = cumCoh/sum(cumCoh(:));
        P = sparse(1:G.N, 1:G.N, 1./sqrt(weight), G.N, G.N);
        
        % Sampling matrix
        ind_obs = datasample(1:G.N, param.nb_meas, ...
            'Replace', true, 'Weights', weight);
        M = sparse(1:param.nb_meas, ind_obs, 1, param.nb_meas, ...
            G.N);
        
        % Generate a signal of unit norm
        x_coeff = randn(param.band_limit, 1);
        x = Uk*x_coeff;
        x = x/norm(x);
        y = M*x;
        e = randn(size(y))*param.noise_list(ind_noise);
        ynoise_init = y + e;
        
        % Prepare matrices for reconstruction
        MP = M*P; MtM = MP'*MP;
        
        % Choose regularisation term
        switch param.reg_term
            case 'L'
                L = G.L;
            case 'L2'
                L = G.L*G.L;
            case 'L4'
                L = G.L*G.L*G.L*G.L;
            case 'Uk_est'
                Uk_est = gsp_estimate_Uk(G, param.band_limit, lk);
        end
        
        % Reconstruct signal for all regularisation parameters
        ynoise = P*P*M'*ynoise_init;
        if strcmp(param.reg_term, 'Uk_est')
            MU = M*Uk_est;
            tic;
            x_est = Uk_est*((MU'*MU)\(MU'*ynoise_init));
            timeRecons(ind_sim, :, ind_noise) = toc;
            err(ind_sim, :, ind_noise) = ...
                norm(x - x_est);
            errUk(ind_sim, :, ind_noise) = ...
                norm(Uk'*(x - x_est));
            errBarUk(ind_sim, :, ind_noise) = ...
                norm((eye(G.N)-Uk*Uk')*(x - x_est));
        else
            for ind_reg = 1:numel(param.reg)
                tic;
                x_est = (MtM + param.reg(ind_reg)*L)\ynoise;
                timeRecons(ind_sim, ind_reg, ind_noise) = toc;
                err(ind_sim, ind_reg, ind_noise) = ...
                    norm(x - x_est);
                errUk(ind_sim, ind_reg, ind_noise) = ...
                    norm(Uk'*(x - x_est));
                errBarUk(ind_sim, ind_reg, ind_noise) = ...
                    norm((eye(G.N)-Uk*Uk')*(x - x_est));
            end
        end
        
    end
end


end

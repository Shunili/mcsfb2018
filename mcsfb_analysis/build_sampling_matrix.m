function [M, selected] = build_sampling_matrix(G, weights, nb_meas, param)
    if nargin < 4
        param = struct;
    end

    if ~isfield(param, 'replacement')
        replacement = 0;
    else
        replacement = param.replacement;
    end

 
    if exist('datasample','builtin')
        if replacement
            ind_obs = datasample(1:G.N, nb_meas, 'Replace', true, 'Weights', weights); 
        else
            ind_obs = datasample(1:G.N, nb_meas, 'Replace', false, 'Weights', weights); 
        end
    elseif exist('RandSampleWR','file') && ~replacement
            if G.N >10000
                warning('This implementation of random sampling without replacement is slow for large graphs. It may be faster to try sampling with replacement.');
            end
            ind_obs = RandSampleWR(1:G.N,nb_meas,weights);
    elseif exist('randsample','file') && replacement
            ind_obs = randsample(1:G.N, nb_meas, true, weights);
    else
        error('Need a function to perform weighted sampling');
    end
    %duplicates=nb_meas-length(unique(ind_obs))
    M = sparse(1:nb_meas, ind_obs, 1, nb_meas, G.N);
    selected=ind_obs';
end

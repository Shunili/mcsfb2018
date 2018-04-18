function [M, selected] = build_sampling_matrix(G, weights, nb_meas, param)
    if nargin < 4
        param = struct;
    end

    if ~isfield(param, 'relacement')
        replacement = 1;
    else
        replacement = param.replacement;
    end

  
    if exist('RandSampleWR')
        if replacement
            ind_obs = randsample(1:G.N, nb_meas, true, weights);
        else
            ind_obs = RandSampleWR(1:G.N,nb_meas,weights);
        end
    elseif exist('datasample')
        if replacement
            ind_obs = datasample(1:G.N, nb_meas, 'Replace', true, 'Weights', weights); 
        else
           ind_obs = datasample(1:G.N, nb_meas, 'Replace', false, 'Weights', weights); 
        end
        
    else
        error('Need a function to perform weighted sampling');
    end
    %duplicates=nb_meas-length(unique(ind_obs))
    M = sparse(1:nb_meas, ind_obs, 1, nb_meas, G.N);
    selected=ind_obs';
end

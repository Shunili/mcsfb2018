function [M, selected] = build_sampling_matrix(G, weights, nb_meas)
    if exist('RandSampleWR')
        ind_obs = RandSampleWR(1:G.N,nb_meas,weights);
    elseif exist('datasample')
        ind_obs = datasample(1:G.N, nb_meas, 'Replace', false, 'Weights', weights); 
    else
        error('Need a function to perform weighted sampling without replacement');
    end
    %duplicates=nb_meas-length(unique(ind_obs))
    M = sparse(1:nb_meas, ind_obs, 1, nb_meas, G.N);
    selected=ind_obs';
end

function [ signal_projections , varargout ] = mcsfb_apply_filters( G, signal, filter_bank, param )

num_bands=size(filter_bank,1);

if ~isfield(param,'order'); param.order = 30; end
if ~isfield(param,'grid_order'); param.grid_order = 1000; end

if ~isfield(param, 'jackson') 
    param.jackson = 0;
end

if ~isfield(param, 'filter_coeffs')

    if param.jackson
        if ~isfield(param, 'shifted_ends')
            error('must pass shifted ends to use jackson cheby coeffs')
        end
        if ~isfield(G, 'lmax')
            G = gsp_estimate_lmax(G);
        end
        param.filter_coeffs = zeros(param.order+1,num_bands);
        for i=1:num_bands
            [~,param.filter_coeffs(:,i)]=gsp_jackson_cheby_coeff(param.shifted_ends(i), param.shifted_ends(i+1),[0 G.lmax], param.order);
        end
    else
        param.filter_coeffs = gsp_cheby_coeff(G, filter_bank, param.order, param.grid_order);
    end 

    varargout{1}=param.filter_coeffs;
end

signal_projections = gsp_vec2mat(gsp_cheby_op(G, param.filter_coeffs, signal, param), num_bands);

end


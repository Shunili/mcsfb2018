function [analysis_coeffs, varargout ]= mcsfb_analysis(G, s, filter_bank, downsampling_sets, param)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

num_bands=size(filter_bank,1);
analysis_coeffs = cell(num_bands, 1);


if gsp_check_fourier(G)
    transformed_coeffs = gsp_vec2mat(gsp_filter_analysis(G,filter_bank,s,param),num_bands);
    varargout{1} = 0; %zeros(param.order+1,num_bands);
else
    if ~isfield(param,'order'); param.order = 30; end
    if ~isfield(param,'grid_order'); param.grid_order = param.order+1; end
    
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
    
    transformed_coeffs = gsp_vec2mat(gsp_cheby_op(G, param.filter_coeffs, s, param), num_bands);
end

for i = 1:num_bands
    analysis_coeffs{i} = transformed_coeffs(downsampling_sets{i},i);
   
end
end
    
% add switch here to use jackson instead of Chebyshev (default in gsp_filter_analysis)?
% if we stay with Chebyshev, add param.grid_pts 

%default for gsp_filter_analysis has order 30, set order & grid_points = 1000
%set param in mcsfb_demo



%add another param in this function: param.jackson = 1. default = 0
%replace gsp_filter_analysis and call gsp_cheby_op(G, coeff, signal)

% In filter bank, store exact filters, add param.approx_coeffs to store the cheby coeffs
% use varargou
% if they are there, call gsp_filter
% elseif jackson is not true, first gsp_cheby_coeff(f) 




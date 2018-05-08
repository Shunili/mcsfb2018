function [analysis_coeffs, varargout ]= mcsfb_analysis(G, s, filter_bank, downsampling_sets, param)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

num_bands=size(filter_bank,1);
analysis_coeffs = cell(num_bands, 1);


if isfield(param,'signal_projections');
    signal_projections=param.signal_projections;
    varargout{1} = 0;
elseif gsp_check_fourier(G)
    signal_projections = gsp_vec2mat(gsp_filter_analysis(G,filter_bank,s,param),num_bands);
    varargout{1} = 0; %zeros(param.order+1,num_bands);
else
    [signal_projections,fc]=mcsfb_apply_filters(G,s,filter_bank,param);
    varargout{1}=fc;
end

for i = 1:num_bands
    analysis_coeffs{i} = signal_projections(downsampling_sets{i},i);
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




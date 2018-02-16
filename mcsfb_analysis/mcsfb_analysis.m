function analysis_coeffs = mcsfb_analysis(G, s, filter_bank, downsampling_sets, param)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

num_bands=size(filter_bank,1);
analysis_coeffs = cell(num_bands, 1);

% add switch here to use jackson instead of Chebyshev (default in gsp_filter_analysis)?
% if we stay with Chebyshev, add param.grid_pts 

transformed_coeffs = gsp_vec2mat(gsp_filter_analysis(G,filter_bank,s,param),num_bands);

for i = 1:num_bands
    analysis_coeffs{i} = transformed_coeffs(downsampling_sets{i},i);
end

end


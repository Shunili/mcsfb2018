function [ f_reconstruct, reconstruction_banded] = mcsfb_sythesis(G, num_bands, downsampling_sets, f_values, shifted_ends, weights, param)
% add one more input argument: param.order
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
if nargin<7
   param = struct;
end

if ~isfield(param,'order')
    param.order=50;
end

reconstruction_banded = cell(num_bands, 1);
f_reconstruct = zeros(G.N, 1);

for i = 1:num_bands
   reconstruction_banded{i} = mcsfb_reconstruct_band(G, downsampling_sets{i}, f_values{i}, shifted_ends(i), shifted_ends(i+1), weights{i}(downsampling_sets{i}));
   f_reconstruct = f_reconstruct + reconstruction_banded{i};
end

end



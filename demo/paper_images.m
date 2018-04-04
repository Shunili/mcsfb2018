close all;
clear all;
randn('seed', 18); 
rand('seed', 18);

% Graph
graph='bunny';

switch graph
    case 'gnp'
        N=100;
        p=.3;
        G=gsp_erdos_renyi(N,p);
    case 'sensor'
        N=500;
        G=gsp_david_sensor_network(N);
    case 'minnesota'
        G=gsp_minnesota(1);
    case 'comet'
        G=gsp_comet(100,20);
    case 'community'
        N=5000;
        G=gsp_community(N);
    case 'bunny'
        G=gsp_bunny();
    otherwise
        error('graph type not recognized');
end

param.compute_full_eigen = 1;


% filter bank
num_bands = 5;
param.band_structure = 0;
param.plot_filters = 1;
param.plot_density_functions = 1;

if ~param.compute_full_eigen
    G = gsp_estimate_lmax(G);
else
    G=gsp_compute_fourier_basis(G);
end

% mcsfb_design_filter_bank
% Design filter bank (pass G, num_bands (M), parameters). Depend on whether
% you have all eigenvalues and what type of filter bank structure you want

%11,01,00,10
param.spacing = 1;
param.spectrum_adapted = 0;

[filter_bank,shifted_ends,band_ends] = mcsfb_design_filter_bank(G,num_bands,param);

% plot filters
if param.plot_filters
    % add eigenvalues to plots?
    figure;
    plot_param.show_sum=0;
    gsp_plot_filter(G,filter_bank,plot_param);
    title('Filters');
end
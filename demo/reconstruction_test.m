G=gsp_bunny();
G=gsp_compute_fourier_basis(G);

%% filter signal here

poly1=1-(G.coords(:,1).^2+G.coords(:,2)-2*G.coords(:,3));
poly2=G.coords(:,1)-G.coords(:,2)+3*G.coords(:,3).^2-2;

signal=zeros(G.N,1);
signal(G.U(:,2)>=0)=poly2(G.U(:,2)>=0);
signal(G.U(:,2)<0)=poly1(G.U(:,2)<0);

tail=(G.A(930,:)==1);
tail(930)=1;
signal(tail)=1;

%% 
param.compute_full_eigen = 0;

if ~param.compute_full_eigen
    G = rmfield(G, 'U');
    G = rmfield(G, 'e');
    G = rmfield(G, 'lmax');
    G = gsp_estimate_lmax(G);
end

% filter bank
num_bands = 4;
param.band_structure = 0;
param.spectrum_adapted=1;
param.plot_filters = 0;
param.plot_density_functions = 0;

% downsampling
param.exact_downsampling_partition=0;
param.plot_downsampling_sets = 0;

% analysis
param.plot_analysis_coeffs = 0;
%param.shifted_ends = shifted_ends;
%param.jackson = 1;

% param.filter_type = 'approximate'; not used right now
% param.filter_type = 'approximate'; not used right now

[filter_bank,shifted_ends,band_ends] = mcsfb_design_filter_bank(G,num_bands,param);

up_limit=shifted_ends(4);
low_limit=shifted_ends(3);
range=[0,G.lmax];
order = 80;
h = @(x) filter_bank{3}(x);
[~, JCH]=gsp_jackson_cheby_coeff(low_limit, up_limit, range, order);
h_tilde = @(x) gsp_cheby_eval(x,JCH,[0,G.lmax]);

f = gsp_filter(G, h_tilde, signal);

%% replacement = 1
%different number of samples compute error, report the average of the errors
param.replacement = 1;

L = ceil(2*log(G.N));

G=spectral_cdf_approx(G, param);
% nb_meas = 20;
% nb_meas=floor((G.spectrum_cdf_approx(up_limit)-G.spectrum_cdf_approx(low_limit))*G.N);

nb_meas = zeros(5,1);
mean_squared_error = zeros(5,1);

for i=1:5
    nb_meas(i) = 200*i;
    [weights, P_min_half] = compute_sampling_weights(G,L,h_tilde);
    [M, selected] = build_sampling_matrix(G, weights, nb_meas(i));
    downsampling_sets = selected;

    %Analysis
    param.order=80;
    transformed_coeffs = gsp_cheby_op(G, JCH, f, param);
    analysis_coeffs = transformed_coeffs(downsampling_sets,1);

    %Sythesis
    f_reconstruct = mcsfb_reconstruct_band2(G, selected, analysis_coeffs, low_limit, up_limit, weights(selected), param);
    error=abs(f-f_reconstruct);
    mean_squared_error(i)=sum(error.^2)/G.N;
end

%% replacement = 0
param.replacement = 0;
L = ceil(2*log(G.N));
G=spectral_cdf_approx(G, param);
% nb_meas = 20;
% nb_meas=floor((G.spectrum_cdf_approx(up_limit)-G.spectrum_cdf_approx(low_limit))*G.N);

nb_meas = zeros(5,1);
mean_squared_error = zeros(5,1);

for i=1:5
    nb_meas(i) = 200*i;
    [weights, P_min_half] = compute_sampling_weights(G,L,h_tilde);
    [M, selected] = build_sampling_matrix(G, weights, nb_meas(i));
    downsampling_sets = selected;

    %Analysis
    param.order=80;
    transformed_coeffs = gsp_cheby_op(G, JCH, f, param);
    analysis_coeffs = transformed_coeffs(downsampling_sets,1);

    %Sythesis
    f_reconstruct = mcsfb_reconstruct_band2(G, selected, analysis_coeffs, low_limit, up_limit, weights(selected), param);
    error=abs(f-f_reconstruct);
    mean_squared_error(i)=sum(error.^2)/G.N;
end

%plot the error here using white black color scheme
%use exact partitioning

%% plot
% plot reconstruction and error for each channel
for i=1:num_bands
    figure;
    param.vertex_size=100;
    param.climits = [-2.5,2.5];
    gsp_plot_signal(G,reconstruction_banded{i}, param);
    caxis([-2.5,2.5]);
    view(0,90)
    title('Reconstruction by Channel');
    set(gca,'FontSize',24);
    
    figure;
    param.vertex_size=100;
    param.climits = [0,2.5];
    gsp_plot_signal(G,abs(gsp_cheby_op(G,filter_coeffs(:,i),f)-reconstruction_banded{i}),param);  
    title('Reconstruction Error by Channel');
    
    view(0,90)
    set(gca,'FontSize',24);
    set(gca,'FontSize',24);
end


%use white to black color scheme
% plot reconstruction
figure;
max_val=max(abs(f));
plot_param.climits = [-max_val,max_val];
plot_param.vertex_size = 100;
gsp_plot_signal(G, f_reconstruct, plot_param);
caxis([-2.5,2.5]);
view(0,90)
set(gca,'FontSize',24);
title('Reconstruction');

% plot reconstruction error
error=abs(f-f_reconstruct);
plot_param.climits = [0, max(error)];
plot_param.vertex_size = 100;
figure;
gsp_plot_signal(G, error, plot_param);
caxis([-2.5,2.5]);
view(0,90);
set(gca,'FontSize',24);
title('Reconstruction Error');

mean_squared_error=sum(error.^2)/G.N
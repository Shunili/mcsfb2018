close all;
clear all;
rand('seed',1);
randn('seed',1);

% Main parameters to explore
num_bands = 4;
order = 50; % used for density estimation and analysis filtering
param.replacement=0;
synth_param.reg_filter=3;
synth_param.order=order; 
synth_param.pcgtol=1e-10;
synth_param.pcgmaxits=2000;
synth_param.gamma=1; % surprisingly insensitive to this parameter % larger gamma puts more weight on matching samples; smaller gamma puts more weight on matching spectral content
param.num_vec=50;
selected_band=3;

% graph and signal
G=gsp_bunny();
load('demo/pwbunny_signal.mat');
signal=pwbunny_signal;
signal=signal-mean(signal);

% Filter bank
G = gsp_estimate_lmax(G);
param.band_structure = 0; % 0: search in the interval
param.spectrum_adapted=1;

[filter_bank, shifted_ends, band_ends,G] = mcsfb_design_filter_bank(G,num_bands,param);

up_limit=shifted_ends(selected_band+1);
low_limit=shifted_ends(selected_band);
range=[0,G.lmax];
h = @(x) filter_bank{selected_band}(x);
[~, JCH]=gsp_jackson_cheby_coeff(low_limit, up_limit, range, order);

h_tilde = @(x) gsp_cheby_eval(x,JCH,[0,G.lmax]);

% Apply polynomial filter to piecewise signal
filtered = gsp_cheby_op(G, JCH, signal);

% Choose number of measurements 
param.cdf_method='kpm';
order=50;
param.order=order;
G=spectral_cdf_approx2(G, param);

%ideal_nb_meas1=round((G.spectrum_cdf_approx(up_limit)-G.spectrum_cdf_approx(low_limit))*G.N)

R=gsp_cheby_opX(G,JCH);
est_num_eigs=gsp_hutch(G,R);
ideal_nb_meas=round(est_num_eigs)
%ideal_nb_meas = 173
%ideal_nb_meas = 445

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orders = 20:10:100;
num_trials=50;
num_m=length(orders);
nmse=zeros(num_m,1);

for i = 1:num_m
    order = orders(i)
    param.order = order; 
    synth_param.order=order;
    total_mse=0;
    %G = rmfield(G,'TkbarLX');
    %G = rmfield(G,'X');
    
    for j = 1:num_trials
        [~, JCH]=gsp_jackson_cheby_coeff(low_limit, up_limit, range, order);
        %h_tilde = @(x) gsp_cheby_eval(x,JCH,[0,G.lmax]);
        
        G = gsp_compute_TkbarLX(G, param);
        R=gsp_cheby_opX(G,JCH);
        norm_Uk= sum(R.^2, 2);
        weights=norm_Uk/sum(norm_Uk);
        [~, selected] = build_sampling_matrix(G, weights, ideal_nb_meas, param);
        analysis_coeffs=filtered(selected(1:ideal_nb_meas));
        [reconstruction, JCH_reg]= mcsfb_reconstruct_band2(G, selected(1:ideal_nb_meas), analysis_coeffs, low_limit, up_limit, weights(selected(1:ideal_nb_meas)), synth_param);
       
        error = abs(reconstruction-filtered);
        mse = sum(error.^2)/G.N;
        nmse(i) = nmse(i) + G.N*mse/sum(filtered.^2);
    end   
end
nmse=nmse/num_trials;


adapted_nmse=zeros(num_m,1);
for i = 1:num_m
    order = orders(i)
    param.order = order; 
    synth_param.order=order;
    total_mse=0;
    %G = rmfield(G,'TkbarLX');
    %G = rmfield(G,'X');
    
    for j = 1:num_trials
        [~, JCH]=gsp_jackson_cheby_coeff(low_limit, up_limit, range, order);
        %h_tilde = @(x) gsp_cheby_eval(x,JCH,[0,G.lmax]);
        
        G = gsp_compute_TkbarLX(G, param);
        R=gsp_cheby_opX(G,JCH);
        norm_Uk= sum(R.^2, 2);
        weights=norm_Uk/sum(norm_Uk);
        adapted_weights=weights.*abs(filtered);
        adapted_weights=adapted_weights/sum(adapted_weights);
        
        [~, selected] = build_sampling_matrix(G, adapted_weights, ideal_nb_meas, param);
        analysis_coeffs=filtered(selected(1:ideal_nb_meas));
        [reconstruction, JCH_reg]= mcsfb_reconstruct_band2(G, selected(1:ideal_nb_meas), analysis_coeffs, low_limit, up_limit, weights(selected(1:ideal_nb_meas)), synth_param);
       
        error = abs(reconstruction-filtered);
        mse = sum(error.^2)/G.N;
        adapted_nmse(i) = adapted_nmse(i) + G.N*mse/sum(filtered.^2);
    end   
end
adapted_nmse=adapted_nmse/num_trials;

figure;
plot1=semilogy(orders,[nmse,adapted_nmse],'-o','LineWidth',3,'MarkerSize',7);
set(plot1, {'MarkerFaceColor'}, get(plot1,'Color')); 
xlabel('Number of Samples');
ylabel('Average Normalized MSE');
set(gca,'FontSize',24);
legend('Not signal adapted','Signal adapted');
grid on;


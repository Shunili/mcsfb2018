close all;
clear all;

rand('seed',0);
randn('seed',0);

num_bands = 4;
%param.order=50; %100; % used for density estimation and analysis filtering
param.search_right_only=0;
param.replacement=0; % use 1 for large graphs for now as the without replacement method is too slow
%param.alpha=.25;
%oversampling_factor=.9; % performance sensitive to this parameter; almost perfect at 1.5
synth_param.reg_filter=3; % reconstruction method: splines
%synth_param.order=param.order; 
synth_param.pcgtol=1e-10; %1e-10; 
synth_param.pcgmaxits=500; %1000; %200; 
synth_param.gamma=1; % surprisingly insensitive to this parameter % larger gamma puts more weight on matching samples; smaller gamma puts more weight on matching spectral content; convergence must faster for larger gamma, which makes sense since we are initializing it to the guesses
%adapted=0; % downsampling sets and number of measurements adapted to signal
extra_plots=0;
param.extra_low_factor=1; % multiplicative factor for extra samples on low channel; taken away from highest channel
param.subtract_mean=0;
param.num_vec=30; % default of 30 for most reconstruction methods seems fine. Plays a more important role if we are trying to reconstruct from the approximated subspaces

% graph and signal
G=gsp_bunny();
load('demo/pwbunny_signal.mat');
signal=pwbunny_signal;
vs=80;
plim=[-2.5,2.5];
initial_signal=signal;
%if param.subtract_mean
%    signal_average=mean(signal);
%    signal=signal-signal_average;
%end

G = gsp_estimate_lmax(G);

param.band_structure = 0; % 0: search in the interval
param.spectrum_adapted=1;


orders=20:10:100;
num_trials=50;
num_m=length(orders);
mean_squared_error = zeros(num_m,1);
adapted_mean_squared_error =  zeros(num_m,1);
adapted = 0;
selected_band = 1;
param.jackson = 1;
            
for k = 1:num_m
    total_mse=0;
        for j = 1:num_trials
            param.order=orders(k);
            synth_param.order=param.order;
            
            [filter_bank, shifted_ends, band_ends, G] = mcsfb_design_filter_bank(G, num_bands, param);
            param.shifted_ends = shifted_ends;
            
            if adapted
                [param.signal_projections,filter_coeffs]=mcsfb_apply_filters(G,signal,filter_bank,param);
                param.adapt_weights=1;
                param.adapt_num_meas=1;  
            end
            
            param.target_samples=G.N-param.subtract_mean;
            [downsampling_sets, weights_banded] = mcsfb_create_downsampling_sets(G, filter_bank, shifted_ends, param);

      
            if adapted % already computed filter coeffs
                [analysis_coeffs] = mcsfb_analysis(G, signal, filter_bank, downsampling_sets, param);
            else
                [analysis_coeffs,filter_coeffs] = mcsfb_analysis(G, signal, filter_bank, downsampling_sets, param);
            end
            
            [f_reconstruct, reconstruction_banded] = mcsfb_synthesis(G, num_bands, downsampling_sets, analysis_coeffs, shifted_ends, weights_banded, synth_param);
            %if param.subtract_mean
            %    f_reconstruct=f_reconstruct+signal_average;
            %end
            
            if ~adapted
                projections=cell(num_bands,1);
                for i=1:num_bands
                    projections{i}=gsp_cheby_op(G,filter_coeffs(:,i),signal);
                end
            else
                projections=mat2cell(param.signal_projections,G.N,ones(1,num_bands));
            end
            
            if adapted
                error=abs(reconstruction_banded{selected_band}-projections(:,selected_band));
            else
                error=abs(reconstruction_banded{selected_band}-projections{selected_band});
            end
            
            total_mse=total_mse+sum(error.^2)/G.N;
        end
    if adapted    
        adapted_mean_squared_error(k) = total_mse/num_trials;
    else
        mean_squared_error(k)=total_mse/num_trials;
    end
end

figure;
hold on;
plot(orders,mean_squared_error,'LineWidth',4);
plot(orders,adapted_mean_squared_error,'LineWidth',4);
xlabel('orders', 'FontSize',24);
ylabel('NMSE', 'FontSize',24);
xlim([20,100]);
set(gca,'FontSize',24);
set(gca,'box','off');
leg = legend('Not Signal Adapted','Signal Adapted');
set(leg,'FontSize',15);
set(leg,'Location','northeast');

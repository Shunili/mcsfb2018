%%
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

% params
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

[filter_bank,shifted_ends,band_ends] = mcsfb_design_filter_bank(G,num_bands,param);

up_limit=shifted_ends(4);
low_limit=shifted_ends(3);
range=[0,G.lmax];
order = 80;
h = @(x) filter_bank{3}(x);
[~, JCH]=gsp_jackson_cheby_coeff(low_limit, up_limit, range, order);
h_tilde = @(x) gsp_cheby_eval(x,JCH,[0,G.lmax]);
f = gsp_cheby_op(G, JCH, signal);

param.vertex_size=100;
figure;
gsp_plot_signal(G,f,param);
view(0,90);

% Check critical sampling value
G=spectral_cdf_approx(G, param);
ideal_nb_meas=floor((G.spectrum_cdf_approx(up_limit)-G.spectrum_cdf_approx(low_limit))*G.N)


%% replacement = 1
%different number of samples compute error, report the average of the errors
param.replacement = 1;

% L = ceil(2*log(G.N));
L=50;

nb_meas=200:100:500;
num_trials=10;
num_m=length(nb_meas);
mean_squared_error = zeros(num_m,1);

[weights, P_min_half] = compute_sampling_weights(G,L,h_tilde);
synth_param.order=order;

for i=1:num_m
    total_mse=0;
    total_error=zeros(G.N,1);
    for j=1:num_trials
        [M, selected] = build_sampling_matrix(G, weights, nb_meas(i),param);

        %Analysis
        analysis_coeffs = f(selected);

        %Sythesis
        f_reconstruct = mcsfb_reconstruct_band2(G, selected, analysis_coeffs, low_limit, up_limit, weights(selected), synth_param);
        error=abs(f-f_reconstruct);
        total_mse=total_mse+sum(error.^2)/G.N;
        total_error=total_error+error;
    end
    mean_squared_error(i)=total_mse/num_trials;
    total_error=total_error/num_trials;
end

figure;
plot(nb_meas,mean_squared_error,'LineWidth',2);

figure;
gsp_plot_signal(G,total_error,param);
view(0,90);
colormap hot;
colormap(flipud(hot))


%% replacement = 0
param.replacement = 0;

[weights, P_min_half] = compute_sampling_weights(G,L,h_tilde);

for i=1:num_m
    total_mse=0;
    total_error=zeros(G.N,1);
    for j=1:num_trials
        [M, selected] = build_sampling_matrix(G, weights, nb_meas(i),param);

        %Analysis
        analysis_coeffs = f(selected);

        %Sythesis
        f_reconstruct = mcsfb_reconstruct_band2(G, selected, analysis_coeffs, low_limit, up_limit, weights(selected), synth_param);
        error=abs(f-f_reconstruct);
        total_mse=total_mse+sum(error.^2)/G.N;
        total_error=total_error+error;
    end
    mean_squared_error(i)=total_mse/num_trials;
    total_error=total_error/num_trials;
end

figure;
plot(nb_meas,mean_squared_error,'LineWidth',2);

figure;
gsp_plot_signal(G,total_error,param);
colormap hot;
colormap(flipud(hot));
view(0,90);


%% reg_filter = 1/(h+e) - 1/(1+e)

param.replacement = 0;
L=50;
nb_meas = ideal_nb_meas;
[weights, P_min_half] = compute_sampling_weights(G,L,h_tilde);
[M, selected] = build_sampling_matrix(G, weights, nb_meas, param);
analysis_coeffs = f(selected);

num_trials=10;
total_mse = 0;
total_error=zeros(G.N,1);
synth_param.reg_filter = 1;
for i=1:num_trials
    f_reconstruct = mcsfb_reconstruct_band2(G, selected, analysis_coeffs, low_limit, up_limit, weights(selected), synth_param);
    error=abs(f-f_reconstruct);
    total_error=total_error+error;
    total_mse=total_mse+sum(error.^2)/G.N;
end

figure;
gsp_plot_signal(G,total_error,param);
colormap hot;
colormap(flipud(hot));
view(0,90);



%% reg_filter = 1-h

% param.replacement = 0;
% L=50;
% nb_meas = ideal_nb_meas;
% [weights, P_min_half] = compute_sampling_weights(G,L,h_tilde);
% [M, selected] = build_sampling_matrix(G, weights, nb_meas, param);
% analysis_coeffs = f(selected);

num_trials=10;
total_mse = 0;
total_error=zeros(G.N,1);
synth_param.reg_filter = 0;
synth_param.precondition = 0;
for i=1:num_trials
    f_reconstruct = mcsfb_reconstruct_band2(G, selected, analysis_coeffs, low_limit, up_limit, weights(selected), synth_param);
    error=abs(f-f_reconstruct);
    total_error=total_error+error;
    total_mse=total_mse+sum(error.^2)/G.N;
end

figure;
gsp_plot_signal(G,total_error/num_trials,param);
colormap hot;
colormap(flipud(hot));
view(0,90);



%% plot reconstruction and error for each channel
% for i=1:num_bands
%     figure;
%     param.climits = [-2.5,2.5];
%     gsp_plot_signal(G,reconstruction_banded{i}, param);
%     caxis([-2.5,2.5]);
%     view(0,90)
%     title('Reconstruction by Channel');
%     set(gca,'FontSize',24);
%     
%     figure;
%     param.vertex_size=100;
%     param.climits = [0,2.5];
%     gsp_plot_signal(G,abs(gsp_cheby_op(G,filter_coeffs(:,i),f)-reconstruction_banded{i}),param);  
%     title('Reconstruction Error by Channel');
%     
%     view(0,90)
%     set(gca,'FontSize',24);
%     set(gca,'FontSize',24);
% end
% 
% 
% %use white to black color scheme
% % plot reconstruction
% figure;
% max_val=max(abs(f));
% plot_param.climits = [-max_val,max_val];
% plot_param.vertex_size = 100;
% gsp_plot_signal(G, f_reconstruct, plot_param);
% caxis([-2.5,2.5]);
% view(0,90)
% set(gca,'FontSize',24);
% title('Reconstruction');
% 
% % plot reconstruction error
% error=abs(f-f_reconstruct);
% plot_param.climits = [0, max(error)];
% plot_param.vertex_size = 100;
% figure;
% gsp_plot_signal(G, error, plot_param);
% caxis([-2.5,2.5]);
% view(0,90);
% set(gca,'FontSize',24);
% title('Reconstruction Error');
% 
% mean_squared_error=sum(error.^2)/G.N
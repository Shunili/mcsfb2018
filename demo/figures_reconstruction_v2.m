%close all;
clear all;
rand('seed',1);
randn('seed',1);
show_spectral=1;

% Main parameters to explore
num_bands = 4;
selected_band=2;
order = 100; % used for density estimation and analysis filtering
oversampling_factor=1; % performance sensitive to this parameter; almost perfect at 1.5
synth_param.reg_filter=3;
synth_param.order=100; 
synth_param.pcgtol=1e-12;
synth_param.pcgmaxits=2000;
synth_param.gamma=1; % surprisingly insensitive to this parameter % larger gamma puts more weight on matching samples; smaller gamma puts more weight on matching spectral content

% Load bunny graph and piecewise smooth signal
G=gsp_bunny();
load('/Users/davidshuman/Dropbox/Current Research Work/MCSFB/Shuni_Thesis/GitHub/mcsfb2018/demo/pwbunny_signal.mat');
signal=pwbunny_signal;

plot_param.vertex_size=100;
plot_param.climits=[-2.5,2.5];
figure;
gsp_plot_signal(G,signal,plot_param);
view(0,90);
set(gca,'FontSize',24);
title('Piecewise smooth signal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% signal definition
% G=gsp_compute_fourier_basis(G);
% poly1=1-(G.coords(:,1).^2+G.coords(:,2)-2*G.coords(:,3));
% poly2=G.coords(:,1)-G.coords(:,2)+3*G.coords(:,3).^2-2;
% 
% signal=zeros(G.N,1);
% signal(G.U(:,2)>=0)=poly2(G.U(:,2)>=0);
% signal(G.U(:,2)<0)=poly1(G.U(:,2)<0);
% 
% tail=(G.A(930,:)==1);
% tail(930)=1;
% signal(tail)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filter bank
G = gsp_estimate_lmax(G);
param.band_structure = 0; 
param.spectrum_adapted=1;

[filter_bank,shifted_ends,band_ends] = mcsfb_design_filter_bank(G,num_bands,param);

% Select and plot one polynomial bandpass filter
up_limit=shifted_ends(selected_band+1);
low_limit=shifted_ends(selected_band);
range=[0,G.lmax];
h = @(x) filter_bank{selected_band}(x);
[~, JCH]=gsp_jackson_cheby_coeff(low_limit, up_limit, range, order);

h_tilde = @(x) gsp_cheby_eval(x,JCH,[0,G.lmax]);
figure;
gsp_plot_filter(G,h_tilde);
set(gca,'FontSize',24);
title('Filter');

% Apply polynomial filter to piecewise signal
filtered = gsp_cheby_op(G, JCH, signal);
plot_lim = max(abs(filtered));
plot_param.climits=[-plot_lim,plot_lim];
% figure;
% gsp_plot_signal(G,filtered,plot_param);
% view(0,90);
% set(gca,'FontSize',24);
% title('Filtered signal');

% Choose number of measurements 
param.cdf_method='kpm';
param.order=order;
param.num_vec=50;
G=spectral_cdf_approx2(G, param);

ideal_nb_meas1=round((G.spectrum_cdf_approx(up_limit)-G.spectrum_cdf_approx(low_limit))*G.N)

R=gsp_cheby_opX(G,JCH);
est_num_eigs=gsp_hutch(G,R);
ideal_nb_meas2=round(est_num_eigs)

% Choose and plot downsampling weights and selected set
norm_Uk= sum(R.^2, 2);
weights=norm_Uk/sum(norm_Uk);
nb_meas=round(oversampling_factor*est_num_eigs)
[~, selected] = build_sampling_matrix(G, weights, nb_meas, param);

plot_param=rmfield(plot_param,'climits');
figure;
gsp_plot_signal(G,weights,plot_param);
colormap(flipud(hot));
view(0,90);
set(gca,'FontSize',24);
title('Non-uniform sampling weights');

sel_verts=zeros(G.N,1);
sel_verts(selected)=1;
figure;
plot_param.climits=[0,2];
gsp_plot_signal(G,sel_verts,plot_param);
colormap(flipud(hot));
view(0,90);
set(gca,'FontSize',24);
title('Selected vertices');

% Analysis coefficients (downsampled on one band)
analysis_coeffs=filtered(selected);
upsampled_analysis_coeffs=zeros(G.N,1);
upsampled_analysis_coeffs(selected)=analysis_coeffs;
plot_param.climits=[-plot_lim,plot_lim];
figure;
gsp_plot_signal(G,upsampled_analysis_coeffs,plot_param);
view(0,90);
set(gca,'FontSize',24);
title('Sample values of filtered signal');

% Reconstruction
[reconstruction, JCH_reg]= mcsfb_reconstruct_band2(G, selected, analysis_coeffs, low_limit, up_limit, weights(selected), synth_param);
figure;
gsp_plot_signal(G,reconstruction,plot_param);
view(0,90); 
set(gca,'FontSize',24);
title('Reconstruction');

plot_param.climits=[-plot_lim,plot_lim];
figure;
gsp_plot_signal(G,filtered,plot_param);
view(0,90);
set(gca,'FontSize',24);
title('Filtered signal');

error=filtered-reconstruction;
mse=sum(error.^2)/G.N
plot_param2.vertex_size=100;
figure;
gsp_plot_signal(G,abs(error),plot_param2);
colormap(flipud(hot));
view(0,90); 
set(gca,'FontSize',24);
title('Reconstruction error (absolute value)');

pen_fun=@(x) gsp_cheby_eval(x,JCH_reg,[0,G.lmax]);
reg_filters = cell(4,1);
reg_filters{1} = @(x) 1-h(x);
reg_filters{2} = @(x) (1-h_tilde(x));
reg_filters{4} = pen_fun;
reg_eps=.01;
reg_filters{3} = @(x) (reg_eps./(h_tilde(x)+reg_eps)-reg_eps/(1+reg_eps));

figure; 
hold on;
ylim([0,2]);
xx = 0:0.001:G.lmax;
plot(xx, h_tilde(xx),'r','LineWidth',3);
plot(xx, reg_filters{1}(xx),'r:','LineWidth',3);
plot(xx, reg_filters{2}(xx),'b','LineWidth',3);
plot(xx, reg_filters{3}(xx),'g','LineWidth',3);
plot(xx, reg_filters{4}(xx),'k','LineWidth',3);
set(gca,'FontSize',24);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);
leg = legend('$h_m(\lambda)$','$1-h_m(\lambda)$','$1-\tilde{h}_m(\lambda)$','$\frac{\epsilon}{\tilde{h}_m(\lambda)+\epsilon}-\frac{\epsilon}{1+\epsilon}$','Damped cubic spline');
set(leg,'Interpreter','latex');
set(leg,'FontSize',20);
set(leg,'Location','southeast');
title('Penalty filters');

% Analyze errors in the spectral domain
if show_spectral
    G2=G;
    G2=gsp_compute_fourier_basis(G2);

    figure;
    gsp_plot_signal_spectral(G2,gsp_gft(G2,filtered),param);
    title('Filtered signal');
    xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);
    ylabel('$|\hat{f}(\lambda)|$','Interpreter','LaTex','FontSize',24);
    set(gca,'FontSize',24);

    figure;
    gsp_plot_signal_spectral(G2,gsp_gft(G2,reconstruction),param);
    title('Reconstruction');
    xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);
    ylabel('$|\hat{f}(\lambda)|$','Interpreter','LaTex','FontSize',24);
    set(gca,'FontSize',24);

    figure;
    gsp_plot_signal_spectral(G2,gsp_gft(G2,error),param);
    title('Reconstruction error');
    xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);
    ylabel('$|\hat{f}(\lambda)|$','Interpreter','LaTex','FontSize',24);
    set(gca,'FontSize',24);
end

% for i = 1:3
%     figure;
%     param.climits = [0 1];
%     gsp_plot_signal_spectral(G,gsp_gft(G,total_error2{i}/num_trials),param);
%     colormap hot;
%     colormap(flipud(hot));
%     view(0,90);
%     xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);
%     ylabel('$|\hat{f}(\lambda)|$','Interpreter','LaTex','FontSize',24);
%     set(gca,'FontSize',24);
% end


% % Repeated trails with different reconstruction methods
% num_trials = 10;
% total_mse = zeros(3,1);
% total_error = cell(3,1);
% total_error2 = cell(3,1);
% for i = 1:3
%     total_mse(i)=0;
%     total_error{i}=zeros(G.N,1);
%     total_error2{i}=zeros(G.N,1);
%     for j = 1:num_trials
%         synth_param.reg_filter = i;
%         synth_param.reg_eps = 1e-2;
%         synth_param.order = 80;
%         reconstruction = mcsfb_reconstruct_band2(G, selected, analysis_coeffs, low_limit, up_limit, weights(selected), synth_param);
%         error=abs(f-reconstruction);
%         error2 = f-reconstruction;
%         total_mse(i)=total_mse(i)+sum(error.^2)/G.N;
%         total_error{i}=total_error{i}+error;
%         total_error2{i}=total_error2{i}+error2;
%     end
% end
% 
% for i = 1:3
%     figure;
%     param.climits = [0 1.5];
%     param.vertex_size=100;
%     gsp_plot_signal(G, total_error{i}/num_trials,param);
%     colormap hot;
%     colormap(flipud(hot));
%     view(0,90);
%     set(gca,'FontSize',24);
%     fig = gcf;
%     fig.InvertHardcopy = 'off';
% end
% 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal-adapted transform

if selected_band>1  % for scaling functions, this could lead to biased interpolation
   % just basing the selection on the coefficients ignores the support of
   % the eigenvectors in the band. We need a combination
   %[~,ii]=sort(abs(filtered),'descend');
   %selected_adapted=ii(1:nb_meas);
   
   weights_adapted=weights.*abs(filtered);
   weights_adapted=weights_adapted/sum(weights_adapted);
   [~, selected_adapted] = build_sampling_matrix(G, weights_adapted, nb_meas, param);
   
   plot_param=rmfield(plot_param,'climits');
   figure;
   gsp_plot_signal(G,weights_adapted,plot_param);
   colormap(flipud(hot));
   view(0,90);
   set(gca,'FontSize',24);
   title('Adapted non-uniform sampling weights');

   sel_verts_adapted=zeros(G.N,1);
   sel_verts_adapted(selected_adapted)=1;
   figure;
   plot_param.climits=[0,2];
   gsp_plot_signal(G,sel_verts_adapted,plot_param);
   view(0,90);
   set(gca,'FontSize',24);
   colormap(flipud(hot));
   title('Adapted selected vertices');
   
   % Reconstruction
   reconstruction_adapted= mcsfb_reconstruct_band2(G, selected_adapted, filtered(selected_adapted), low_limit, up_limit, weights_adapted(selected_adapted), synth_param);
   plot_param.climits=[-plot_lim,plot_lim];  
   figure;
   gsp_plot_signal(G,reconstruction_adapted,plot_param);
   view(0,90); 
   set(gca,'FontSize',24);
   title('Adapted reconstruction');

   error_adapted=filtered-reconstruction_adapted;
   mse_adapted=sum(error_adapted.^2)/G.N
   plot_param2.vertex_size=100;
   figure;
   gsp_plot_signal(G,abs(error_adapted),plot_param2);
   colormap(flipud(hot));
   view(0,90); 
   set(gca,'FontSize',24);
   title('Adapted reconstruction error (absolute value)');

   % Analyze adapted errors in the spectral domain
   figure;
   gsp_plot_signal_spectral(G2,gsp_gft(G2,reconstruction_adapted),param);
   title('Adapted reconstruction');
   xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);
   ylabel('$|\hat{f}(\lambda)|$','Interpreter','LaTex','FontSize',24);
   set(gca,'FontSize',24);

   figure;
   gsp_plot_signal_spectral(G2,gsp_gft(G2,error_adapted),param);
   title('Adapted reconstruction error');
   xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);
   ylabel('$|\hat{f}(\lambda)|$','Interpreter','LaTex','FontSize',24);
   set(gca,'FontSize',24);
end




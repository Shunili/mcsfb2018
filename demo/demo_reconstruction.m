close all;
clear all;
rand('seed',1);
randn('seed',1);
show_spectral=1;
extra_plots=0;
for_paper=0;
extra_sim=0;

% Main parameters to explore
num_bands = 4;
selected_band=1;
order = 100; % used for density estimation and analysis filtering
param.replacement=0;
synth_param.reg_filter=3;
synth_param.order=order; 
synth_param.pcgtol=1e-10;
synth_param.pcgmaxits=2000;
synth_param.gamma=1; % surprisingly insensitive to this parameter % larger gamma puts more weight on matching samples; smaller gamma puts more weight on matching spectral content
param.num_vec=50;


% Load bunny graph and piecewise smooth signal
G=gsp_bunny();
load('demo/pwbunny_signal.mat');
signal=pwbunny_signal;
signal=signal-mean(signal);

plot_param.vertex_size=100;
if extra_plots
    plot_param.climits=[-2.5,2.5];
    figure;
    gsp_plot_signal(G,signal,plot_param);
    view(0,90);
    set(gca,'FontSize',24);
    if ~for_paper
        title('Piecewise smooth signal');
    end
end

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

[filter_bank,shifted_ends,band_ends,G] = mcsfb_design_filter_bank(G,num_bands,param);

% Select and plot one polynomial bandpass filter
up_limit=shifted_ends(selected_band+1);
low_limit=shifted_ends(selected_band);
range=[0,G.lmax];
h = @(x) filter_bank{selected_band}(x);
[~, JCH]=gsp_jackson_cheby_coeff(low_limit, up_limit, range, order);

h_tilde = @(x) gsp_cheby_eval(x,JCH,[0,G.lmax]);

if extra_plots
    figure;
    gsp_plot_filter(G,h_tilde);
    xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);
    set(gca,'FontSize',24);
    title('Filter');
end

% Apply polynomial filter to piecewise signal
filtered = gsp_cheby_op(G, JCH, signal);
plot_lim = max(abs(filtered));
plot_param.climits=[-plot_lim,plot_lim];
figure;
gsp_plot_signal(G,filtered,plot_param);
view(0,90);
set(gca,'FontSize',24);
if ~for_paper
    title('Filtered signal');
end

% Choose number of measurements 
param.cdf_method='kpm';
param.order=order;
G=spectral_cdf_approx2(G, param);

%ideal_nb_meas1=round((G.spectrum_cdf_approx(up_limit)-G.spectrum_cdf_approx(low_limit))*G.N)

R=gsp_cheby_opX(G,JCH);
est_num_eigs=gsp_hutch(G,R);
ideal_nb_meas=round(est_num_eigs)

% Choose and plot downsampling weights and selected set
norm_Uk= sum(R.^2, 2);
weights=norm_Uk/sum(norm_Uk);

plot_param=rmfield(plot_param,'climits');
figure;
gsp_plot_signal(G,weights,plot_param);
colormap(flipud(hot));
view(0,90);
set(gca,'FontSize',24);
if ~for_paper
    title('Non-uniform sampling weights');
end

if extra_sim
    num_factors=41; 
    num_trials=50; 
else 
    num_factors=9;
    num_trials=5;
end

oversampling_factors=linspace(1,3,num_factors);
nb_meas=round(oversampling_factors*est_num_eigs)
nmse=zeros(num_factors,1);
error=zeros(G.N,num_factors,num_trials);

for k=1:num_trials
    [~, selected] = build_sampling_matrix(G, weights, nb_meas(end), param);
    if mod(k,10)==0
        mess=sprintf('Trial %d',k);
        display(mess);
    end
    for i=1:num_factors
        analysis_coeffs=filtered(selected(1:nb_meas(i)));
        [reconstruction, JCH_reg]= mcsfb_reconstruct_band2(G, selected(1:nb_meas(i)), analysis_coeffs, low_limit, up_limit, weights(selected(1:nb_meas(i))), synth_param);
        error(:,i,k)=filtered-reconstruction;
        mse=sum(error(:,i,k).^2)/G.N;
        nmse(i)=nmse(i)+G.N*mse/sum(filtered.^2);
    end
end
nmse=nmse/num_trials;

adapted_nmse=zeros(num_factors,1);
adapted_error=zeros(G.N,num_factors,num_trials);
adapted_weights=weights.*abs(filtered);
adapted_weights=adapted_weights/sum(adapted_weights);

figure;
gsp_plot_signal(G,adapted_weights,plot_param);
colormap(flipud(hot));
view(0,90);
set(gca,'FontSize',24);
if ~for_paper
    title('Signal-adapted non-uniform sampling weights');
end

for k=1:num_trials
    [~, adapted_selected] = build_sampling_matrix(G, adapted_weights, nb_meas(end), param);
    if mod(k,10)==0
        mess=sprintf('Adapted Trial %d',k);
        display(mess);
    end
    for i=1:num_factors
        adapted_analysis_coeffs=filtered(adapted_selected(1:nb_meas(i)));
        [adapted_reconstruction, JCH_reg]= mcsfb_reconstruct_band2(G, adapted_selected(1:nb_meas(i)), adapted_analysis_coeffs, low_limit, up_limit, adapted_weights(adapted_selected(1:nb_meas(i))), synth_param);
        adapted_error(:,i,k)=filtered-adapted_reconstruction;
        adapted_mse=sum(adapted_error(:,i,k).^2)/G.N;
        adapted_nmse(i)=adapted_nmse(i)+G.N*adapted_mse/sum(filtered.^2);
    end
end
adapted_nmse=adapted_nmse/num_trials;

figure;
plot1=semilogy(nb_meas,[nmse,adapted_nmse],'-o','LineWidth',3,'MarkerSize',7);
set(plot1, {'MarkerFaceColor'}, get(plot1,'Color')); 
xlabel('Number of Samples');
ylabel('Average Normalized MSE');
set(gca,'FontSize',24);
legend('Not signal adapted','Signal adapted');
grid on;

selected_factor=1;
%if ~for_paper
    sel_verts=zeros(G.N,1);
    sel_verts(selected(1:nb_meas(selected_factor)))=1;
    figure;
    plot_param.climits=[0,2];
    gsp_plot_signal(G,sel_verts,plot_param);
    colormap(flipud(hot));
    view(0,90);
    set(gca,'FontSize',24);
    colorbar off;
    if ~for_paper
        title('Selected vertices');
    end
    
    adapted_sel_verts=zeros(G.N,1);
    adapted_sel_verts(adapted_selected(1:nb_meas(selected_factor)))=1;
    figure;
    plot_param.climits=[0,2];
    gsp_plot_signal(G,adapted_sel_verts,plot_param);
    colormap(flipud(hot));
    view(0,90);
    set(gca,'FontSize',24);
    if ~for_paper
        title('Selected vertices (signal-adapted)');
    end
    colorbar off;
%end

if extra_plots
    % Analysis coefficients (downsampled on one band)
    upsampled_analysis_coeffs=zeros(G.N,1);
    upsampled_analysis_coeffs(selected)=analysis_coeffs;
    plot_param.climits=[-plot_lim,plot_lim];
    figure;
    gsp_plot_signal(G,upsampled_analysis_coeffs,plot_param);
    view(0,90);
    set(gca,'FontSize',24);
    title('Sample values of filtered signal');

    % Reconstruction
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
end

selected_trial=1;
plot_param2.vertex_size=100;
if extra_plots
    emax=max(abs([error(:,selected_factor,selected_trial);adapted_error(:,selected_factor,selected_trial)]));
    plot_param2.climits=[0,emax];
    figure;
    gsp_plot_signal(G,abs(error(:,selected_factor,selected_trial)),plot_param2);
    colormap(flipud(hot));
    view(0,90); 
    set(gca,'FontSize',24);
    if ~for_paper
        title('Reconstruction error (absolute value)');
    end

    figure;
    gsp_plot_signal(G,abs(adapted_error(:,selected_factor,selected_trial)),plot_param2);
    colormap(flipud(hot));
    view(0,90); 
    set(gca,'FontSize',24);
    if ~for_paper
        title('Signal-adapted reconstruction error (absolute value)');
    end
    
    aemax0=max(abs([mean((squeeze(error(:,selected_factor,:))),2);mean((squeeze(adapted_error(:,selected_factor,:))),2)]));
    plot_param2.climits=[-aemax0,aemax0];
    figure;
    gsp_plot_signal(G,mean((squeeze(error(:,selected_factor,:))),2),plot_param2);
    view(0,90); 
    set(gca,'FontSize',24);
    if ~for_paper
        title('Average reconstruction error');
    end
    
    figure;
    gsp_plot_signal(G,mean((squeeze(adapted_error(:,selected_factor,:))),2),plot_param2);
    view(0,90); 
    set(gca,'FontSize',24);
    if ~for_paper
        title('Average signal-adapted reconstruction error');
    end
end

avg_error=mean(abs(squeeze(error(:,selected_factor,:))),2);
avg_error_adapted=mean(abs(squeeze(adapted_error(:,selected_factor,:))),2);
aemax=max([avg_error;avg_error_adapted]);
plot_param2.climits=[0,aemax];
figure;
gsp_plot_signal(G,avg_error,plot_param2);
colormap(flipud(hot));
view(0,90); 
set(gca,'FontSize',24);
if ~for_paper
    title('Average reconstruction error (absolute value)');
end

figure;
gsp_plot_signal(G,avg_error_adapted,plot_param2);
colormap(flipud(hot));
view(0,90); 
set(gca,'FontSize',24);
if ~for_paper
    title('Average signal-adapted reconstruction error (absolute value)');
end




pen_fun=@(x) gsp_cheby_eval(x,JCH_reg,[0,G.lmax]);
reg_filters = cell(4,1);
reg_filters{1} = @(x) 1-h(x);
reg_filters{2} = @(x) (1-h_tilde(x));
reg_filters{4} = pen_fun;
reg_eps=(sqrt(5)-1)/2;
reg_filters{3} = @(x) (1./(h_tilde(x)+reg_eps)-1/(1+reg_eps));

if ~for_paper
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
    leg = legend('$h_m(\lambda)$','$1-h_m(\lambda)$','$1-\tilde{h}_m(\lambda)$','$\frac{1}{\tilde{h}_m(\lambda)+\epsilon}-\frac{1}{1+\epsilon}$','Damped cubic spline');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',20);
    set(leg,'Location','southeast');
    title('Penalty filters');
end

% Analyze errors in the spectral domain
if show_spectral
    G2=G;
    G2=gsp_compute_fourier_basis(G2);

    figure;
    gsp_plot_signal_spectral(G2,gsp_gft(G2,filtered),param);
    if ~for_paper
        title('Filtered signal');
    end
    xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);
    ylabel('$|\hat{f}(\lambda)|$','Interpreter','LaTex','FontSize',24);
    set(gca,'FontSize',24);

    if extra_plots
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

% if selected_band>1  % for scaling functions, this could lead to biased interpolation
%    % just basing the selection on the coefficients ignores the support of
%    % the eigenvectors in the band. We need a combination
%    %[~,ii]=sort(abs(filtered),'descend');
%    %selected_adapted=ii(1:nb_meas);
%    
%    weights_adapted=weights.*abs(filtered);
%    weights_adapted=weights_adapted/sum(weights_adapted);
%    [~, selected_adapted] = build_sampling_matrix(G, weights_adapted, nb_meas, param);
%    
%    plot_param=rmfield(plot_param,'climits');
%    figure;
%    gsp_plot_signal(G,weights_adapted,plot_param);
%    colormap(flipud(hot));
%    view(0,90);
%    set(gca,'FontSize',24);
%    title('Adapted non-uniform sampling weights');
% 
%    sel_verts_adapted=zeros(G.N,1);
%    sel_verts_adapted(selected_adapted)=1;
%    figure;
%    plot_param.climits=[0,2];
%    gsp_plot_signal(G,sel_verts_adapted,plot_param);
%    view(0,90);
%    set(gca,'FontSize',24);
%    colormap(flipud(hot));
%    title('Adapted selected vertices');
%    
%    % Reconstruction
%    reconstruction_adapted= mcsfb_reconstruct_band2(G, selected_adapted, filtered(selected_adapted), low_limit, up_limit, weights_adapted(selected_adapted), synth_param);
%    plot_param.climits=[-plot_lim,plot_lim];  
%    figure;
%    gsp_plot_signal(G,reconstruction_adapted,plot_param);
%    view(0,90); 
%    set(gca,'FontSize',24);
%    title('Adapted reconstruction');
% 
%    error_adapted=filtered-reconstruction_adapted;
%    mse_adapted=sum(error_adapted.^2)/G.N
%    plot_param2.vertex_size=100;
%    figure;
%    gsp_plot_signal(G,abs(error_adapted),plot_param2);
%    colormap(flipud(hot));
%    view(0,90); 
%    set(gca,'FontSize',24);
%    title('Adapted reconstruction error (absolute value)');
% 
%    % Analyze adapted errors in the spectral domain
%    figure;
%    gsp_plot_signal_spectral(G2,gsp_gft(G2,reconstruction_adapted),param);
%    title('Adapted reconstruction');
%    xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);
%    ylabel('$|\hat{f}(\lambda)|$','Interpreter','LaTex','FontSize',24);
%    set(gca,'FontSize',24);
% 
%    figure;
%    gsp_plot_signal_spectral(G2,gsp_gft(G2,error_adapted),param);
%    title('Adapted reconstruction error');
%    xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);
%    ylabel('$|\hat{f}(\lambda)|$','Interpreter','LaTex','FontSize',24);
%    set(gca,'FontSize',24);
% end
% 
% 
% 

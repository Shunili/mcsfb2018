G=gsp_bunny();
G=gsp_compute_fourier_basis(G);

poly1=1-(G.coords(:,1).^2+G.coords(:,2)-2*G.coords(:,3));
poly2=G.coords(:,1)-G.coords(:,2)+3*G.coords(:,3).^2-2;

signal=zeros(G.N,1);
signal(G.U(:,2)>=0)=poly2(G.U(:,2)>=0);
signal(G.U(:,2)<0)=poly1(G.U(:,2)<0);

tail=(G.A(930,:)==1);
tail(930)=1;
signal(tail)=1;

ee = G.e;
U = G.U;

% params
param.compute_full_eigen = 0;

if ~param.compute_full_eigen
    G = rmfield(G, 'U');
    G = rmfield(G, 'e');
    G = rmfield(G, 'lmax');
    G = gsp_estimate_lmax(G);
end

num_bands = 4;
param.band_structure = 0; 
param.spectrum_adapted=1;

[filter_bank,shifted_ends,band_ends] = mcsfb_design_filter_bank(G,num_bands,param);

up_limit=shifted_ends(4);
low_limit=shifted_ends(3);
range=[0,G.lmax];
order = 80;
h = @(x) filter_bank{3}(x);
[~, JCH]=gsp_jackson_cheby_coeff(low_limit, up_limit, range, order);
h_tilde = @(x) gsp_cheby_eval(x,JCH,[0,G.lmax]);
f = gsp_cheby_op(G, JCH, signal);

%%%%
param.vertex_size=100;
param.climits=[-2,2];
figure;
gsp_plot_signal(G,f,param);
view(0,90);
set(gca,'FontSize',24);

G.e = ee;
G.U = U;
figure;
gsp_plot_signal_spectral(G,gsp_gft(G,f),param);
view(0,90);
% title('signal');
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);
ylabel('$|\hat{f}(\lambda)|$','Interpreter','LaTex','FontSize',24);
set(gca,'FontSize',24);


% ideal_nb_meas=floor((G.spectrum_cdf_approx(3.5)-G.spectrum_cdf_approx(1.2))*G.N)
G=spectral_cdf_approx(G, param);
ideal_nb_meas=floor((G.spectrum_cdf_approx(up_limit)-G.spectrum_cdf_approx(low_limit))*G.N)

%%%%

param.replacement = 0;
L=50;
nb_meas = ideal_nb_meas;
[weights, P_min_half] = compute_sampling_weights(G,L,h_tilde);
[M, selected] = build_sampling_matrix(G, weights, nb_meas, param);
analysis_coeffs = f(selected);

G = rmfield(G, 'U');
G = rmfield(G, 'e');
G = rmfield(G, 'lmax');
G = gsp_estimate_lmax(G);

num_trials = 10;
total_mse = zeros(3,1);
total_error = cell(3,1);
total_error2 = cell(3,1);
for i = 1:3
    total_mse(i)=0;
    total_error{i}=zeros(G.N,1);
    total_error2{i}=zeros(G.N,1);
    for j = 1:num_trials
        synth_param.reg_filter = i;
        synth_param.reg_eps = 1e-2;
        synth_param.order = 80;
        f_reconstruct = mcsfb_reconstruct_band2(G, selected, analysis_coeffs, low_limit, up_limit, weights(selected), synth_param);
        error=abs(f-f_reconstruct);
        error2 = f-f_reconstruct;
        total_mse(i)=total_mse(i)+sum(error.^2)/G.N;
        total_error{i}=total_error{i}+error;
        total_error2{i}=total_error2{i}+error2;
    end
end

for i = 1:3
    figure;
    param.climits = [0 1.5];
    param.vertex_size=100;
    gsp_plot_signal(G, total_error{i}/num_trials,param);
    colormap hot;
    colormap(flipud(hot));
    view(0,90);
    set(gca,'FontSize',24);
    fig = gcf;
    fig.InvertHardcopy = 'off';
end

G.e = ee;
G.U = U;

for i = 1:3
    figure;
    param.climits = [0 1];
    gsp_plot_signal_spectral(G,gsp_gft(G,total_error2{i}/num_trials),param);
    colormap hot;
    colormap(flipud(hot));
    view(0,90);
    xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);
    ylabel('$|\hat{f}(\lambda)|$','Interpreter','LaTex','FontSize',24);
    set(gca,'FontSize',24);
end
G = rmfield(G, 'U');
G = rmfield(G, 'e');
G = rmfield(G, 'lmax');
G = gsp_estimate_lmax(G);
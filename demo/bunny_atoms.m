close all;
clear all;
rand('seed',1);
randn('seed',1);


% Main parameters to explore
num_bands = 5;
selected_band=3;
param.order = 30; % used for density estimation and analysis filtering
param.replacement=0;
param.num_vec=50;
param.subtract_mean=0;
adapted=0;



% Load bunny graph and piecewise smooth signal
G=gsp_bunny();
load('/Users/davidshuman/Dropbox/Current_Research_Work/MCSFB/Shuni_Thesis/GitHub/mcsfb2018/demo/pwbunny_signal.mat');
signal=pwbunny_signal;
signal=signal-mean(signal);

% Filter bank
G = gsp_estimate_lmax(G);
param.band_structure = 0; 
param.spectrum_adapted=1;

[filter_bank,shifted_ends,band_ends,G] = mcsfb_design_filter_bank(G,num_bands,param);
[~,filter_coeffs]=gsp_jackson_cheby_coeff(shifted_ends(selected_band),shifted_ends(selected_band+1),[0,G.lmax],param.order);

figure;
xx=0:.01:G.lmax;
plot(xx,gsp_cheby_eval(xx,filter_coeffs,[0,G.lmax]),'r','LineWidth',3);
set(gca,'FontSize',24);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);
ylabel('$\tilde{h}_m(\lambda)$','Interpreter','LaTex','FontSize',24);

% create downsampling sets (and ensure critical sampling)
param.jackson = 1;
param.shifted_ends = shifted_ends;

if adapted
    [param.signal_projections,filter_coeffs]=mcsfb_apply_filters(G,signal,filter_bank,param);
    param.adapt_weights=1;
    param.adapt_num_meas=1;  
end

[downsampling_sets, weights_banded] = mcsfb_create_downsampling_sets(G, filter_bank, shifted_ends, param);

% atoms at selected band
num_atoms=length(downsampling_sets{selected_band});
X=sparse(downsampling_sets{selected_band},1:num_atoms,ones(num_atoms,1),G.N,num_atoms);

atoms=gsp_cheby_op(G,filter_coeffs,X);
%atoms=gsp_filter(G,filter_bank{selected_band},X,param);

% plot atom
atom_ind=1;
figure;
param.vertex_size=100;
maxval=max(abs(atoms(:,atom_ind)));
param.climits=[-maxval,maxval];
gsp_plot_signal(G,atoms(:,atom_ind),param);
set(gca,'FontSize',24);
view(0,90);

% plot all atoms in spectral domain
G=gsp_compute_fourier_basis(G);
atoms_hat=G.U'*atoms;

figure;
plot(G.e,abs(atoms_hat));
hold on;
plot(G.e, mean(abs(atoms_hat),2),'k','LineWidth',3);
set(gca,'FontSize',24);
box on;
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);



figure;
stem(G.e,abs(atoms_hat(:,atom_ind)));
set(gca,'FontSize',24);
box on;
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);
ylabel('$|\tilde{h}_m({\cal L})U^{\top}\delta_i|$','Interpreter','LaTex','FontSize',24);



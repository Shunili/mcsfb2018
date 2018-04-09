close all;
clear all;
randn('seed', 18); 
rand('seed', 18);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Sensor Network
close all;
clear all;
randn('seed', 18); 
rand('seed', 18);

G=gsp_david_sensor_network(500);
param.compute_full_eigen = 1;
% filter bank
num_bands = 3;
param.band_structure = 0;
param.plot_filters = 1;
param.plot_density_functions = 1;

if ~param.compute_full_eigen
    G = gsp_estimate_lmax(G);
else
    G=gsp_compute_fourier_basis(G);
end

%11,01,00,10
param.spacing = 1;
param.spectrum_adapted = 1;

[filter_bank,shifted_ends,band_ends] = mcsfb_design_filter_bank(G,num_bands,param);
% 

figure;
plot_param.plot_eigenvalues = 1;
plot_param.show_sum=0;
gsp_plot_filter(G,filter_bank, plot_param);
set(gca,'FontSize',24);
xlim([0,G.lmax+1]);
set(gca,'box','off');
% set(gca, 'XTick', [0,G.e',15]);
% set(gca,'XTickLabel',[]);
set(gca, 'YTick', 0:0.5:1);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);

order = 200;
filter_coeffs = zeros(order+1,num_bands);
for i=1:num_bands
    [~,filter_coeffs(:,i)]=gsp_jackson_cheby_coeff(shifted_ends(i), shifted_ends(i+1),[0 G.lmax], order);
end

approx_filters=cell(num_bands,1);
for i=1:num_bands
    approx_filters{i}=@(x) gsp_cheby_eval(x,filter_coeffs(:,i),[0,G.lmax]);
end

figure;
plot_param.show_sum=0;
plot_param.plot_eigenvalues = 1;
gsp_plot_filter(G,approx_filters,plot_param);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24) 
set(gca,'FontSize',24);
xlim([0,G.lmax+1]);
set(gca,'box','off');
% set(gca, 'XTick', [0,G.e',15]);
% set(gca,'xticklabel',{[]}) 
set(gca, 'YTick', 0:0.5:1);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% New Graph
load('net25.mat');
A=Problem.A;
A = A - diag(diag(A)); 
G=gsp_graph(A);

G=gsp_compute_fourier_basis(G);

shifted_ends(1) = 0;
shifted_ends(2) = 10;
shifted_ends(3) = 40;
shifted_ends(4) = 166.8325;

for l = 1:3
    filter_bank{l}=@(x) ((shifted_ends(l) <= x) & (x <= shifted_ends(l+1)));
end


figure;
plot_param.plot_eigenvalues = 1;
plot_param.show_sum=0;
gsp_plot_filter(G,filter_bank, plot_param);
set(gca,'FontSize',24);
xlim([0,G.lmax+1]);
set(gca,'box','off');
% set(gca, 'XTick', [0,C',15]);
% set(gca,'XTickLabel',[]);
set(gca, 'YTick', 0:0.5:1);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);

order = 200;
filter_coeffs = zeros(order+1,num_bands);
for i=1:num_bands
    [~,filter_coeffs(:,i)]=gsp_jackson_cheby_coeff(shifted_ends(i), shifted_ends(i+1),[0 G.lmax], order);
end

approx_filters=cell(num_bands,1);
for i=1:num_bands
    approx_filters{i}=@(x) gsp_cheby_eval(x,filter_coeffs(:,i),[0,G.lmax]);
end

figure;
plot_param.show_sum=0;
plot_param.plot_eigenvalues = 1;
gsp_plot_filter(G,approx_filters,plot_param);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24) 
set(gca,'FontSize',24);
xlim([0,G.lmax+1]);
set(gca,'box','off');
% set(gca, 'XTick', [0,G.e',15]);
% set(gca,'xticklabel',{[]}) 
set(gca, 'YTick', 0:0.5:1);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);

%%%%%%%%%%%%%%%%%%%%%%%%

param.order = 30;
if ~isfield(G,'spectrum_cdf_approx')
    [G, cdf_vals]= spectral_cdf_approx(G, param);
end

xx = 0:0.001:G.lmax;
yy_cdf = G.spectrum_cdf_approx(xx);

param.compute_full_eigen = 0;
G = gsp_estimate_lmax(G);

% filter bank
num_bands = 3;
param.band_structure = 0;
param.plot_filters = 0;
param.spectrum_adapted=1;
param.plot_density_func1ions = 0;

% downsampling
param.exact_downsampling_partition=0;
param.plot_downsampling_sets = 0;

% analysis
param.plot_analysis_coeffs = 0;
%param.shifted_ends = shifted_ends;
%param.jackson = 1;


[filter_bank, shifted_ends, band_ends] = mcsfb_design_filter_bank(G, num_bands,param);


G=gsp_compute_fourier_basis(G);


figure;
plot(xx, yy_cdf,'b','LineWidth',4);
hold on;
scatter(band_ends, [0,1/4,1/2,1],'r','LineWidth',4);
grid on;
box on;
xlim([0,G.lmax+1]);
ylim([0,1.02]);
xlabel('\lambda');
set(gca,'XTick',0:30:(G.lmax+1));
set(gca,'FontSize',24);


order = 200;
filter_coeffs = zeros(order+1,num_bands);
for i=1:num_bands
    [~,filter_coeffs(:,i)]=gsp_jackson_cheby_coeff(band_ends(i), band_ends(i+1),[0 G.lmax], order);
end

approx_filters=cell(num_bands,1);
for i=1:num_bands
    approx_filters{i}=@(x) gsp_cheby_eval(x,filter_coeffs(:,i),[0,G.lmax]);
end

figure;
plot_param.show_sum=0;
plot_param.plot_eigenvalues = 1;
gsp_plot_filter(G,approx_filters,plot_param);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24) 
set(gca,'FontSize',24);
xlim([0,G.lmax+1]);
set(gca,'box','off');
% set(gca, 'XTick', [0,G.e',15]);
% set(gca,'xticklabel',{[]}) 
set(gca, 'YTick', 0:0.5:1);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);


figure;
plot(xx, yy_cdf,'b','LineWidth',4);
hold on;
scatter(band_ends, [0,1/4,1/2,1],'r','LineWidth',4);
scatter(shifted_ends, [0,1/4,1/2,1],'g','LineWidth',4);
grid on;
box on;
xlim([0,G.lmax+1]);
ylim([0,1.02]);
xlabel('\lambda');
set(gca,'XTick',0:30:(G.lmax+1));
set(gca,'FontSize',24);


order = 200;
filter_coeffs = zeros(order+1,num_bands);
for i=1:num_bands
    [~,filter_coeffs(:,i)]=gsp_jackson_cheby_coeff(shifted_ends(i), shifted_ends(i+1),[0 G.lmax], order);
end

approx_filters=cell(num_bands,1);
for i=1:num_bands
    approx_filters{i}=@(x) gsp_cheby_eval(x,filter_coeffs(:,i),[0,G.lmax]);
end

figure;
plot_param.show_sum=0;
plot_param.plot_eigenvalues = 1;
gsp_plot_filter(G,approx_filters,plot_param);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24) 
set(gca,'FontSize',24);
xlim([0,G.lmax+1]);
set(gca,'box','off');
% set(gca, 'XTick', [0,G.e',15]);
% set(gca,'xticklabel',{[]}) 
set(gca, 'YTick', 0:0.5:1);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);







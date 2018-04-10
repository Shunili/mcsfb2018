%% Sensor Network
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

param.spacing = 1;
param.spectrum_adapted = 1;
[filter_bank,shifted_ends,band_ends] = mcsfb_design_filter_bank(G,num_bands,param);

order = 80;
filter_coeffs = zeros(order+1,num_bands);
for i=1:num_bands
    [~,filter_coeffs(:,i)]=gsp_jackson_cheby_coeff(shifted_ends(i), shifted_ends(i+1),[0 G.lmax], order);
end
approx_filters=cell(num_bands,1);
for i=1:num_bands
    approx_filters{i}=@(x) gsp_cheby_eval(x,filter_coeffs(:,i),[0,G.lmax]);
end

hh = cell(2,1);
hh{1} = filter_bank{2};
hh{2} = approx_filters{2};
figure;
plot_param.show_sum=0;
plot_param.plot_eigenvalues = 1;
plot_param.x_tic = 2;
gsp_plot_filter(G,hh,plot_param);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24) 
set(gca,'FontSize',24);
xlim([0,G.lmax+1]);
set(gca,'box','off');
% set(gca, 'XTick', [0,G.e',15]);
% set(gca,'xticklabel',{[]}) 
set(gca, 'YTick', 0:0.5:1);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Net25 Graph
load('net25.mat');
A=Problem.A;
A = A - diag(diag(A)); 
A(4228,6327) = 1;
A(6327,4228) = 1;
G=gsp_graph(A);

shifted_ends(1) = 0;
shifted_ends(2) = 10;
shifted_ends(3) = 40;
shifted_ends(4) = 166.8325;

for l = 1:3
    filter_bank{l}=@(x) ((shifted_ends(l) <= x) & (x <= shifted_ends(l+1)));
end

order = 80;
filter_coeffs = zeros(order+1,num_bands);
for i=1:num_bands
    [~,filter_coeffs(:,i)]=gsp_jackson_cheby_coeff(shifted_ends(i), shifted_ends(i+1),[0 G.lmax], order);
end
approx_filters=cell(num_bands,1);
for i=1:num_bands
    approx_filters{i}=@(x) gsp_cheby_eval(x,filter_coeffs(:,i),[0,G.lmax]);
end

hh = cell(2,1);
hh{1} = filter_bank{2};
hh{2} = approx_filters{2};

figure;
plot_param.show_sum=0;
plot_param.plot_eigenvalues = 1;
plot_param.x_tic = 30;
gsp_plot_filter(G,hh,plot_param);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24) 
set(gca,'FontSize',24);
xlim([0,170]);
set(gca,'box','off');
set(gca, 'YTick', 0:0.5:1);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);

%% Filter Bank Design (Net25 demo)

param.order = 30;
if ~isfield(G,'spectrum_cdf_approx')
    [G, cdf_vals]= spectral_cdf_approx(G, param);
end

if ~isfield(G,'spectrum_pdf_approx')
    xx = 0:0.001:G.lmax;
    delta=.1;
    G.spectrum_pdf_approx = @(x) (G.spectrum_cdf_approx(x+delta) - G.spectrum_cdf_approx(x-delta)) / (2*delta);% first derivative
end

xx = 0:0.001:G.lmax;
yy_pdf = G.spectrum_pdf_approx(xx);
yy_cdf = G.spectrum_cdf_approx(xx);

param.compute_full_eigen = 0;
if ~param.compute_full_eigen
    G = gsp_estimate_lmax(G);
else
    G=gsp_compute_fourier_basis(G);
end

% filter bank
num_bands = 3;
param.band_structure = 0;
param.plot_filters = 0;
param.spectrum_adapted=1;
param.plot_density_func1ions = 0;

[filter_bank, shifted_ends, band_ends] = mcsfb_design_filter_bank(G, num_bands,param);

G=gsp_compute_fourier_basis(G); % for plotting x ticks

figure;
plot(xx, yy_cdf,'k','LineWidth',3);
hold on;
scatter(band_ends, [0,1/4,1/2,1],100,'b','filled');
a=annotation('line',[0.13,.26],[.535,.535]);
a.Color='red';
a.LineStyle=':';
a.LineWidth=2;
b=annotation('line',[.275,.275],[.53,.16]);
b.Color='red';
b.LineStyle=':';
b.LineWidth=2;
c=annotation('line',[0.13,.21],[.35,.35]);
c.Color='red';
c.LineStyle=':';
c.LineWidth=2;
d=annotation('line',[.225,.225],[.345,.16]);
d.Color='red';
d.LineStyle=':';
d.LineWidth=2;
grid on;
box on;
xlim([0,170]);
ylim([0,1.02]);
xlabel('\lambda');
set(gca,'XTick',0:30:(G.lmax+1));
set(gca,'YTick',0:0.25:1);
set(gca,'FontSize',24);

order = 80;
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
plot_param.x_tic = 30;
gsp_plot_filter(G,approx_filters,plot_param);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24) 
set(gca,'FontSize',24);
xlim([0, 170]);
set(gca,'box','off');
% set(gca, 'XTick', [0,unique(G.e)',170], 'sz', 50);
% set(gca,'xticklabel',{}) 
set(gca, 'YTick', 0:0.5:1);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);

L = zeros(4,1);
L(1) = 0;
L(4) = 0;
U = zeros(4,1);
U(1) = 0;
U(4) = 0;
  for k = 2:3
    L(k) = band_ends(k)-(band_ends(k)+band_ends(k-1))/2;
    U(k) = (band_ends(k)+band_ends(k+1))/2-band_ends(k);
  end
  
figure;
plot(xx, yy_pdf,'k','LineWidth',3);
hold on;
errorbar(band_ends, G.spectrum_pdf_approx(band_ends),L,U,'horizontal', 'bx','MarkerSize',16,'LineWidth',4);
scatter(shifted_ends, G.spectrum_pdf_approx(shifted_ends),100,'r','filled');
a=annotation('arrow',[.30,.35],[.64,.22]);
a.Color='red';
a.LineStyle=':';
a.LineWidth=2;
b=annotation('arrow',[.24,.21],[.47,.25]);
b.Color='red';
b.LineStyle=':';
b.LineWidth=2;
grid on;
box on;
xlim([0,170]);
ylim([0,0.035]);
xlabel('\lambda');
set(gca,'XTick',0:30:170);
set(gca,'FontSize',24);

order = 80;
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
plot_param.x_tick = 30;
gsp_plot_filter(G,approx_filters,plot_param);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24) 
set(gca,'FontSize',24);
xlim([0,G.lmax+1]);
set(gca,'box','off');
set(gca, 'YTick', 0:0.5:1);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);







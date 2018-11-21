clear all;
close all;

rand('seed',0);
randn('seed',0);

% Initialize graph
graph_type='sensor500';
switch graph_type
    case 'minnesota'
        G=gsp_minnesota(1);
    case 'sensor64'
        G=gsp_sensor(64);
    case 'sensor500'
        G=gsp_david_sensor_network(500);
    otherwise
        error('Unknown graph type');
end

G = gsp_compute_fourier_basis(G);
 
% create and plot signal f here 
poly1a=-2*G.coords(:,1)+.5;
poly2a=G.coords(:,1).^2+G.coords(:,2).^2+.5;
f=zeros(G.N,1);
p1=(G.coords(:,2)>=(1-G.coords(:,1))) & (G.coords(:,2)<(1.5-G.coords(:,1)));
p2=(G.coords(:,2)<(0.6-G.coords(:,1)));
p3= p1 | p2;
f(~p3)=poly1a(~p3);
f(p3)=poly2a(p3);

plot_param.vertex_size=100;
figure;
gsp_plot_signal(G,f,plot_param);
title('Original Signal');

% Design filter bank
num_bands = 5;
param.band_structure = 0;
param.plot_filters = 1;
[filter_bank,shifted_ends] = mcsfb_design_filter_bank(G,num_bands,param);

xx=0:.001:ceil(G.lmax);
filter_data=zeros(length(xx),num_bands);
for i=1:num_bands
    filter_data(:,i)=filter_bank{i}(xx);
end

figure;
set(gca, 'ColorOrder', [1 0 0;1 0.5 .2; 1 1 0; 0 1 0; 0 0 1], 'NextPlot', 'replacechildren');
plot(xx,filter_data,'LineWidth',3);
set(gca,'FontSize',24);
ax = gca;
ax.YTick = [0 0.5 1];
set(gca, 'XTick', [0,G.e',15]);
xTickLabels = cell(1,G.N+2);  % Empty cell array the same length as xAxis
xTickLabels{1} = 0;
xTickLabels{G.N+2}=15;
                                 % Fills in only the values you want
set(gca,'XTickLabel',xTickLabels);   % Update the tick labels
box off;
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24) 

% figure;
% plot_param.show_sum=0;
% gsp_plot_filter(G,filter_bank,plot_param);

% Identify and plot partition into uniqueness sets
param.exact_downsampling_partition=1;
param.reverse=0;
[downsampling_sets, weights_banded] = mcsfb_create_downsampling_sets(G, filter_bank, shifted_ends, param);

partitionIDs=zeros(G.N,1);
for i=1:num_bands
    partitionIDs(downsampling_sets{i})=i;
end

figure;
param.vertex_size=100;
gsp_plot_signal(G,partitionIDs,param);
hcb=colorbar;
colormap(prism(5));
caxis([.5,5.5]);
set(hcb,'YTick',1:5);
set(hcb,'FontSize',24);


% Test if we found uniqueness sets
total_samples=0;
chosen=zeros(G.N,1);
cond_nums=zeros(num_bands,1);
for i=1:num_bands
    total_samples=total_samples+length(downsampling_sets{i});
    chosen(downsampling_sets{i})=1;
    cond_nums(i)=cond(G.U(downsampling_sets{i},filter_bank{i}(G.e)>0));
end
cond_nums=cond_nums

if (total_samples==G.N && sum(chosen)==G.N)
    display('Uniqueness partition successfully created');
end

% mcsfb_analysis
analysis_coeffs = mcsfb_analysis(G, f, filter_bank, downsampling_sets, param);

% plot coeffs by channel
param.plot_analysis_coeffs = 1;
if param.plot_analysis_coeffs
    for i=1:num_bands
        figure;
        chan_coeffs=zeros(G.N,1);
        chan_coeffs(downsampling_sets{i})=analysis_coeffs{i};
        gsp_plot_signal(G,chan_coeffs,plot_param);
        caxis([-max(abs(chan_coeffs)),max(abs(chan_coeffs))]);
    end
end

% mcsfb_sythesis
[f_reconstruct, reconstruction_banded] = mcsfb_synthesis(G, num_bands, downsampling_sets, analysis_coeffs, shifted_ends, weights_banded,param);

% plot reconstruction
figure;
max_val=max(abs(f));
plot_param.climits = [-max_val,max_val];
gsp_plot_signal(G, f_reconstruct, plot_param);
title('Reconstruction');

% plot reconstruction error
error=abs(f-f_reconstruct);
plot_param.climits = [0, max(error)];
figure;
gsp_plot_signal(G, error, plot_param);
title('Reconstruction Error');

% plot reconstruction error by band
for i=1:num_bands
    figure;
    gsp_plot_signal(G, abs(reconstruction_banded{i}-gsp_filter(G,filter_bank{i},f)), plot_param);
    title('Reconstruction Error by Band');
end


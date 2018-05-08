close all;
clear all;

rand('seed',0);
randn('seed',0);


MyData=csvread('/Users/davidshuman/Dropbox/Current Research Work/MCSFB/Shuni_Thesis/GitHub/mcsfb2018/Data/MyData2.csv');
avg_temp=MyData(:,1);
coords=MyData(:,2:3);
inds=MyData(:,4);

% create 8 neighbor graph and find subset where we have data
param.eight_neighbor=1;
G0=gsp_2dgrid(1385,596,param);
W=G0.W(inds,inds);

% remove disconnected nodes
disconnected=(sum(W,2)==0);
inds=inds(~disconnected);
avg_temp=avg_temp(~disconnected);
coords=coords(~disconnected,:);
W1=G0.W(inds,inds);

% remove small components
G1=gsp_graph(W1,coords);
[V,D]=eigs(G1.L+1e-11*speye(G1.N),4,'sm'); % make sure this works for other machines. compute total number in each and choose the large one
large_component=(abs(V(:,4))>.000000001);
W=G1.W(large_component,large_component);
coords=coords(large_component,:);
avg_temp=avg_temp(large_component);

% create graph and plot signal
G=gsp_graph(W,coords);
is_connected=gsp_check_connectivity(G)

figure;
param.vertex_size=1;
param.climits=[-12,25];
gsp_plot_signal(G,avg_temp,param);
set(gca,'FontSize',24);


% design filter bank
G = gsp_estimate_lmax(G);
num_bands = 6;
param.band_structure = 0;
param.spectrum_adapted=1;
param.order=80;

tic
[filter_bank,shifted_ends,band_ends,G] = mcsfb_design_filter_bank(G,num_bands,param);
filter_design_time=toc

% plot cdf
figure;
xxx=0:.001:G.lmax;
plot(xxx, G.spectrum_cdf_approx(xxx),'k','LineWidth',3);
hold on;
plot(xxx,xxx/G.lmax,'k:','LineWidth',1.8);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24) 
set(gca,'FontSize',24);
ylim([0,1]);
xlim([0,G.lmax]);
set(gca,'XTick',0:2:G.lmax);

% plot filters
figure;
plot_param.show_sum=0;
gsp_plot_filter(G,filter_bank,plot_param);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24) 
set(gca,'FontSize',24);

% create downsampling sets (and ensure critical sampling)
adapted=1;
if adapted
    param.signal_projections=mcsfb_apply_filters(G,signal,filter_bank,param);
    param.adapt_weights=0;
end

tic
[downsampling_sets, weights_banded] = mcsfb_create_downsampling_sets(G, filter_bank, shifted_ends, param);


% ensure critical sampling
total_samples=0;
for i=1:num_bands
    total_samples=total_samples+length(downsampling_sets{i});
end
if total_samples>G.N % eliminate from last band
    extra=total_samples-G.N;
    new_length=length(downsampling_sets{num_bands})-extra;
    downsampling_sets{num_bands}=downsampling_sets{num_bands}(1:new_length);
elseif total_samples<G.N % resample first band with more samples
    additional=G.N-total_samples;
    num_first=length(downsampling_sets{1})+additional;
    [~, selected] = build_sampling_matrix(G, weights_banded{1}, num_first);
    downsampling_sets{1}=selected;
end
downsampling_selection_time=toc

% Plot sampling weights for each band
sampling_param.climits=[0,max((cell2mat(weights_banded)))];
sampling_param.vertex_size=70;
for i=1:num_bands
    figure;
    gsp_plot_signal(G,(weights_banded{i}),sampling_param);
    %titlestr=sprintf('Sampling weights for channel %d',i);
    %title(titlestr);
    colormap(flipud(hot));
    set(gca,'FontSize',24);
    set(gcf,'color',[211,211,211]/255);
end

% plot samples chosen for first band
figure;
band1_param.vertex_size=0.5;
band1=zeros(G.N,1);
band1(downsampling_sets{1})=1;
gsp_plot_signal(G,band1,band1_param);
colormap(flipud(hot));
colorbar off;

% perform analysis
param.shifted_ends = shifted_ends;
param.jackson = 1;

tic
[analysis_coeffs,filter_coeffs] = mcsfb_analysis(G, avg_temp, filter_bank, downsampling_sets, param);
analysis_time=toc


% plot all analysis coeffs
all_analysis_coeff_mag=[abs(analysis_coeffs{1}),ones(size(analysis_coeffs{1}));
    abs(analysis_coeffs{2}),2*ones(size(analysis_coeffs{2}));
    abs(analysis_coeffs{3}),3*ones(size(analysis_coeffs{3}));
    abs(analysis_coeffs{4}),4*ones(size(analysis_coeffs{4}));
    abs(analysis_coeffs{5}),5*ones(size(analysis_coeffs{5}));
    abs(analysis_coeffs{6}),6*ones(size(analysis_coeffs{6}))];
[~,ii]=sort(all_analysis_coeff_mag(:,1),'descend');
all_analysis_coeff_mag_sorted=all_analysis_coeff_mag(ii,:);

figure;
scatter(1:G.N,all_analysis_coeff_mag_sorted(:,1),5,all_analysis_coeff_mag_sorted(:,2));
box on;
set(gca,'FontSize',24);

% plot wavelet coeffs
all_wavelet_coeff_mag=sort(abs(cell2mat(analysis_coeffs(2:num_bands))),'descend');
figure;
scatter(1:length(all_wavelet_coeff_mag),all_wavelet_coeff_mag);
box on;
set(gca,'FontSize',24);


% Plot approx filters
approx_filters=cell(num_bands,1);
for i=1:num_bands
   approx_filters{i}=@(x) gsp_cheby_eval(x,filter_coeffs(:,i),[0,G.lmax]);
end
figure;
gsp_plot_filter(G,approx_filters,plot_param); 
title('Approximate Filters');
set(gca,'FontSize',24);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);

% plot filtered signal by channel
projections=cell(num_bands,1);
for i=1:num_bands
    projections{i}=gsp_cheby_op(G,filter_coeffs(:,i),avg_temp);
    figure;
    gsp_plot_signal(G,projections{i},param);
    %titlestr=sprintf('Channel %d Filtered Signal',i);
    %title(titlestr);
    set(gca,'FontSize',24);
end

proj2=gsp_cheby_op(G,filter_coeffs(:,2),avg_temp);
param.climits=[0,max(abs(proj2))];
figure;
gsp_plot_signal(G,abs(proj2),param);
title('Absolute value of channel 2 filtered signal');
colormap(flipud(hot));

% plot downsampled coeffs by channel
for i=1:num_bands
    figure;
    chan_coeffs=zeros(G.N,1);
    chan_coeffs(downsampling_sets{i})=analysis_coeffs{i};
    plotlim=max(abs(chan_coeffs));
    param.climits=[-plotlim,plotlim];
    gsp_plot_signal(G,chan_coeffs,param);
    titlestr=sprintf('Channel %d Downsampled Analysis Coefficients',i);
    title(titlestr);
end


% plot an example atom from each band
atom_param.vertex_size=1;
for i=1:num_bands
    figure;
    atom=gsp_cheby_op(G,filter_coeffs(:,i),gsp_delta(G,downsampling_sets{i}(1)));
    atom_param.climits=[-max(abs(atom)),max(abs(atom))];
    gsp_plot_signal(G,atom,atom_param);
    %titlestr=sprintf('Atom from Channel %d',i);
    %title(titlestr);
    set(gca,'FontSize',24);
end

% % salt lake atoms
% slc_inds=find((G.coords(:,1)>-112.5) & (G.coords(:,1)<-112.3) & (G.coords(:,2)<41.5) & (G.coords(:,2)>41));
% slc_band1=intersect(slc_inds,downsampling_sets{1});
% 
% 
% atom=gsp_cheby_op(G,filter_coeffs(:,1),gsp_delta(G,slc_band1(1)));
% atom_param.climits=[-max(abs(atom)),max(abs(atom))];
% figure;
% gsp_plot_signal(G,atom,atom_param);
% set(gca,'FontSize',24);



% plot zoom-in of Mass
ma_verts=((G.coords(:,1)>-71.7)&(G.coords(:,1)<-69)&(G.coords(:,2)>36)&(G.coords(:,2)<43));
sum(ma_verts)
WMA=G.W(ma_verts,ma_verts);
coordsMA=G.coords(ma_verts,:);
GMA=gsp_graph(WMA,coordsMA);
GMA.plotting.vertex_size=20;
GMA.plotting.edge_color='k';
figure;
gsp_plot_graph(GMA);

bos_verts=((G.coords(:,1)>-71.4)&(G.coords(:,1)<-70.5)&(G.coords(:,2)>42)&(G.coords(:,2)<42.8));
sum(bos_verts)
WBos=G.W(bos_verts,bos_verts);
coordsBos=G.coords(bos_verts,:);
GBos=gsp_graph(WBos,coordsBos);
GBos.plotting.vertex_size=20;
GBos.plotting.edge_color='k';
figure;
gsp_plot_graph(GBos);

% synthesis
synth_param.order=100;
[f_reconstruct, reconstruction_banded] = mcsfb_sythesis(G, num_bands, downsampling_sets, analysis_coeffs, shifted_ends, weights_banded, synth_param);

% plot overall reconstruction
param.climits=[-12,25];
figure;
gsp_plot_signal(G,f_reconstruct,param);
set(gca,'FontSize',24);

% plot overall reconstruction error
error_param.vertex_size=1;
figure;
gsp_plot_signal(G,abs(f_reconstruct-avg_temp),error_param);
set(gca,'FontSize',24);

% plot reconstruction by band
for i=1:num_bands
    figure;
    gsp_plot_signal(G,reconstruction_banded{i},param);
    set(gca,'FontSize',24);
end

% plot reconstruction error by band
error_param.vertex_size=1;
for i=1:num_bands
    figure;
    gsp_plot_signal(G,abs(reconstruction_banded{i}-projections{i}),error_param);
    set(gca,'FontSize',24);
end

close all;
clear all;

rand('seed',0);
randn('seed',0);

%TODO: change number of samples, update reconsruct_band2 (new method,speed
%up), add spectral plots

% Main parameters to explore
num_bands = 4;
param.order=100; % used for density estimation and analysis filtering
param.search_right_only=0;
param.replacement=0; % use 1 for large graphs for now as the without replacement method is too slow
%oversampling_factor=.9; % performance sensitive to this parameter; almost perfect at 1.5
synth_param.reg_filter=3; % reconstruction method: splines
synth_param.order=param.order; 
synth_param.pcgtol=1e-10; 
synth_param.pcgmaxits=500; 
synth_param.gamma=1; % surprisingly insensitive to this parameter % larger gamma puts more weight on matching samples; smaller gamma puts more weight on matching spectral content; convergence must faster for larger gamma, which makes sense since we are initializing it to the guesses
adapted=1; % downsampling sets adapted to signal
extra_plots=0;
param.extra_low_factor=2; % multiplicative factor for extra samples on low channel; taken away from highest channel
param.subtract_mean=1;
param.num_vec=30; % default of 30 for most reconstruction methods seems fine. Plays a more important role if we are trying to reconstruct from the approximated subspaces

% Graph and signal
graph='bunny';

switch graph
    case 'gnp'
        N=100;
        p=.3;
        G=gsp_erdos_renyi(N,p);
    case 'sensor'
        N=500;
        G=gsp_david_sensor_network(N);
    case 'minnesota'
        G=gsp_minnesota(1);
        signal=(G.coords(:,1)>-93);
        vs=80;
        plim=[-1.5,1.5];
    case 'comet'
        G=gsp_comet(100,20);
    case 'bunny'
        G=gsp_bunny();
        load('/Users/davidshuman/Dropbox/Current Research Work/MCSFB/Shuni_Thesis/GitHub/mcsfb2018/demo/pwbunny_signal.mat');
        signal=pwbunny_signal;
        vs=80;
        plim=[-2.5,2.5];
    case 'community'
        N=5000;
        G=gsp_community(N);
    case 'net25'
        load('/Users/davidshuman/Dropbox/Current Research Work/MCSFB/Shuni_Thesis/GitHub/mcsfb2018/net25 graph data/net25.mat');
        A=Problem.A;
        A = A - diag(diag(A)); 
        A(4228,6327) = 1;
        A(6327,4228) = 1;
        G=gsp_graph(A);
    case 'temperature'
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

        % create graph and signal
        G=gsp_graph(W,coords);
        signal=avg_temp;
        vs=1;
        plim=[-12,25];
    otherwise
        error('graph type not recognized');
end

is_connected=gsp_check_connectivity(G)

figure;
param.vertex_size=vs;
param.climits=plim;
gsp_plot_signal(G,signal,param);
set(gca,'FontSize',24);
if strcmp(graph,'bunny')
    view(0,90);
end
title('Initial signal');

initial_signal=signal;
if param.subtract_mean
    signal_average=mean(signal)
    signal=signal-signal_average;
    figure;
    param.vertex_size=vs;
    param.climits=plim;
    gsp_plot_signal(G,signal,param);
    set(gca,'FontSize',24);
    if strcmp(graph,'bunny')
        view(0,90);
    end
    title('Initial mean normalized signal');
end

% design filter bank
G = gsp_estimate_lmax(G);
param.band_structure = 0;
param.spectrum_adapted=1;

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
title('Estimated spectral cumulative distribution function');

% plot filters
plot_param.show_sum=0;
if extra_plots
    figure;
    gsp_plot_filter(G,filter_bank,plot_param);
    xlabel('$\lambda$','Interpreter','LaTex','FontSize',24) 
    set(gca,'FontSize',24);
end

% create downsampling sets (and ensure critical sampling)
param.jackson = 1;
param.shifted_ends = shifted_ends;
tic
if adapted
    [param.signal_projections,filter_coeffs]=mcsfb_apply_filters(G,signal,filter_bank,param);
    param.adapt_weights=1;
end
filter_time=toc

tic
[downsampling_sets, weights_banded] = mcsfb_create_downsampling_sets(G, filter_bank, shifted_ends, param);


% ensure critical sampling
total_samples=0;
target_samples=G.N-param.subtract_mean;
for i=1:num_bands
    total_samples=total_samples+length(downsampling_sets{i});
end
if total_samples>target_samples % eliminate from last band
    extra=total_samples-target_samples;
    new_length=length(downsampling_sets{num_bands})-extra;
    downsampling_sets{num_bands}=downsampling_sets{num_bands}(1:new_length);
elseif total_samples<target_samples % resample first band with more samples
    additional=target_samples-total_samples;
    num_first=length(downsampling_sets{1})+additional;
    [~, selected] = build_sampling_matrix(G, weights_banded{1}, num_first);
    downsampling_sets{1}=selected;
end
downsampling_selection_time=toc

% Plot sampling weights for each band
sampling_param.climits=[0,max((cell2mat(weights_banded)))];
sampling_param.vertex_size=vs;
for i=1:num_bands
    figure;
    gsp_plot_signal(G,(weights_banded{i}),sampling_param);
    titlestr=sprintf('Sampling weights for channel %d',i);
    title(titlestr);
    colormap(flipud(hot));
    set(gca,'FontSize',24);
    set(gcf,'color',[211,211,211]/255);
    if strcmp(graph,'bunny')
        view(0,90);
    end
end

% plot samples chosen for first band
figure;
band1_param.vertex_size=vs/2;
band1_param.climits=[0,2];
band1=zeros(G.N,1);
band1(downsampling_sets{1})=1;
gsp_plot_signal(G,band1,band1_param);
colorbar off;
title('Selected samples for channel 1');
colormap(flipud(hot));
if strcmp(graph,'bunny')
    view(0,90);
end
set(gca,'FontSize',24);


% perform analysis
tic
if adapted % already computed filter coeffs
    [analysis_coeffs] = mcsfb_analysis(G, signal, filter_bank, downsampling_sets, param);
else
    [analysis_coeffs,filter_coeffs] = mcsfb_analysis(G, signal, filter_bank, downsampling_sets, param);
end
analysis_time=toc

if extra_plots
    % plot all analysis coeffs
    all_analysis_coeff_mag=[];
    for i=1:num_bands
        all_analysis_coeff_mag=[all_analysis_coeff_mag;[abs(analysis_coeffs{i}),i*ones(size(analysis_coeffs{i}))]];
    end
    [~,ii]=sort(all_analysis_coeff_mag(:,1),'descend');
    all_analysis_coeff_mag_sorted=all_analysis_coeff_mag(ii,:);

    figure;
    scatter(1:G.N,all_analysis_coeff_mag_sorted(:,1),5,all_analysis_coeff_mag_sorted(:,2));
    box on;
    set(gca,'FontSize',24);
    title('All coeffs');

    % plot wavelet coeffs
    all_wavelet_coeff_mag=sort(abs(cell2mat(analysis_coeffs(2:num_bands))),'descend');
    figure;
    scatter(1:length(all_wavelet_coeff_mag),all_wavelet_coeff_mag);
    box on;
    set(gca,'FontSize',24);
    title('Wavelet coeffs');
end

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
if ~adapted
    projections=cell(num_bands,1);
    for i=1:num_bands
        projections{i}=gsp_cheby_op(G,filter_coeffs(:,i),signal);
    end
else
    projections=mat2cell(param.signal_projections,G.N,ones(1,num_bands));
end

param.climits=plim;
for i=1:num_bands
    figure;
    if i>1
        maxval=max(abs(projections{i}));
        param.climits=[-maxval,maxval];
    end
    gsp_plot_signal(G,projections{i},param);
    titlestr=sprintf('Channel %d Filtered Signal',i);
    title(titlestr);
    set(gca,'FontSize',24);
    if strcmp(graph,'bunny')
        view(0,90);
    end
end

if extra_plots
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
        set(gca,'FontSize',24);
        if strcmp(graph,'bunny')
            view(0,90);
        end
    end


    % plot an example atom from each band
    atom_param.vertex_size=vs;
    for i=1:num_bands
        figure;
        atom=gsp_cheby_op(G,filter_coeffs(:,i),gsp_delta(G,downsampling_sets{i}(1)));
        atom_param.climits=[-max(abs(atom)),max(abs(atom))];
        gsp_plot_signal(G,atom,atom_param);
        titlestr=sprintf('Atom from Channel %d',i);
        title(titlestr);
        set(gca,'FontSize',24);
        if strcmp(graph,'bunny')
            view(0,90);
        end
    end
end


% synthesis
tic
[f_reconstruct, reconstruction_banded] = mcsfb_sythesis(G, num_bands, downsampling_sets, analysis_coeffs, shifted_ends, weights_banded, synth_param);
if param.subtract_mean
    f_reconstruct=f_reconstruct+signal_average;
end
synth_time=toc

% plot overall reconstruction

param.climits=plim;
figure;
gsp_plot_signal(G,f_reconstruct,param);
set(gca,'FontSize',24);
title('Reconstruction');
if strcmp(graph,'bunny')
    view(0,90);
end

% plot overall reconstruction error
error=f_reconstruct-initial_signal;
mse=sum(error.^2)/G.N
error_param.vertex_size=vs;
figure;
gsp_plot_signal(G,abs(error),error_param);
set(gca,'FontSize',24);
title('Reconstruction error');
if strcmp(graph,'bunny')
    view(0,90);
end
colormap(flipud(hot));

% plot reconstruction by band
param.climits=plim;
for i=1:num_bands
    if i>1
        maxval=max(abs(projections{i}));
        param.climits=[-maxval,maxval];
    end
    figure;
    gsp_plot_signal(G,reconstruction_banded{i},param);
    set(gca,'FontSize',24);
    titlestr=sprintf('Reconstruction from Channel %d',i);
    title(titlestr);
    if strcmp(graph,'bunny')
        view(0,90);
    end
end

% plot reconstruction error by band
mse_by_band=zeros(num_bands,1);
for i=1:num_bands
    figure;
    gsp_plot_signal(G,abs(reconstruction_banded{i}-projections{i}),error_param);
    set(gca,'FontSize',24);
    titlestr=sprintf('Reconstruction error for Channel %d',i);
    title(titlestr);
    if strcmp(graph,'bunny')
        view(0,90);
    end
    colormap(flipud(hot));
    mse_by_band(i)=sum((reconstruction_banded{i}-projections{i}).^2)/G.N;
end
mse_by_band

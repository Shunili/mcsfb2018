clear all;
close all;

rand('seed',0);
randn('seed',0);

% Initialize graph
G=gsp_bunny();
G=gsp_compute_fourier_basis(G);

% create and plot signal f here 
% create signal
poly1=1-(G.coords(:,1).^2+G.coords(:,2)-2*G.coords(:,3));
poly2=G.coords(:,1)-G.coords(:,2)+3*G.coords(:,3).^2-2;

signal=zeros(G.N,1);
signal(G.U(:,2)>=0)=poly2(G.U(:,2)>=0);
signal(G.U(:,2)<0)=poly1(G.U(:,2)<0);

tail=(G.A(930,:)==1);
tail(930)=1;
signal(tail)=1;
f = signal;

% plot signal in both domains
param.vertex_size=100;
figure;
gsp_plot_signal(G,signal,param);
caxis([-2.5,2.5]);
view(0,90)
title('Signal in the Vertex Domain');

figure;
gsp_plot_signal_spectral(G,gsp_gft(G,signal));
set(gca,'FontSize',24);
xlabel('$\lambda$','Interpreter','latex','FontSize',24);
ylabel('$|\hat{f}(\lambda)|$','Interpreter','latex','FontSize',24);
title('Signal in the Spectral Domain');

% Params

% general
param.compute_full_eigen = 0;

% filter bank
num_bands = 5;
param.band_structure = 0;
param.spectrum_adapted=1;
param.plot_filters = 1;
param.plot_density_functions = 1;

% downsampling
param.exact_downsampling_partition=0;
param.plot_downsampling_sets = 0;

% analysis
param.plot_analysis_coeffs = 0;
%param.shifted_ends = shifted_ends;
%param.jackson = 1;

% param.filter_type = 'approximate'; not used right now
if ~param.compute_full_eigen
    G = rmfield(G, 'U');
    G = rmfield(G, 'e');
    G = rmfield(G, 'lmax');
    G = gsp_estimate_lmax(G);
end

% mcsfb_design_filter_bank
% Design filter bank (pass G, num_bands (M), parameters). Depend on whether
% you have all eigenvalues and what type of filter bank structure you want

[filter_bank,shifted_ends,band_ends] = mcsfb_design_filter_bank(G,num_bands,param);

% plot filters
if param.plot_filters
    figure;
    plot_param.show_sum=0;
    gsp_plot_filter(G,filter_bank,plot_param);
    title('Filters');
end

% plot cdf, pdf, eigenvalues, band_ends
if param.plot_density_functions
    if ~isfield(G,'spectrum_cdf_approx')
        [G, cdf_vals]= spectral_cdf_approx(G, param);
    end

    if ~isfield(G,'spectrum_pdf_approx')
      xx = 0:0.001:G.lmax;
      delta=.001;
      G.spectrum_pdf_approx = @(x) (G.spectrum_cdf_approx(x+delta) - G.spectrum_cdf_approx(x-delta)) / (2*delta);% first derivative
    end
    xx = 0:0.001:G.lmax;
    yy_cdf = G.spectrum_cdf_approx(xx);
    yy_pdf = G.spectrum_pdf_approx(xx);
%     
%     if ~isfield(G,'e')
%         G = gsp_compute_fourier_basis(G);
%     end
    figure; 
    hold on;
    plot(xx, yy_cdf,'LineWidth',3);
    plot(xx, yy_pdf,'LineWidth',3);
    ylim([0,1]);
    set(gca,'FontSize',24);
    ax = gca;
    ax.YTick = [0 0.5 1];
%     ax.XTick = [0,uniqG.e',15];
%     ax.XTickLabel = [];
    box off;
    xlabel('$\lambda$','Interpreter','LaTex','FontSize',24) 
end

% mcsfb_create_downsampling_sets
% Choose downsampling sets (either randomly or deterministically, depending
% on whether you have U; pass G and filter cell to this function). Output
% is a cell of length M where each entry has a vector with downsample
% vertices

[downsampling_sets, weights_banded] = mcsfb_create_downsampling_sets(G, filter_bank, shifted_ends, param);

% plot selected vertices

if param.plot_downsampling_sets
    if param.exact_downsampling_partition
        partition_ids=zeros(G.N,1);
        for i=1:num_bands
            partition_ids(downsampling_sets{i})=i;
        end
        plot_param.climits=[1,num_bands];
        figure;
        gsp_plot_signal(G,partition_ids,plot_param);
        hcb=colorbar;
        colormap(prism(num_bands));
        caxis([.5,num_bands+.5]);
        set(hcb,'YTick',1:5);
        set(hcb,'FontSize',24);
        title('Partition');
    else
        for i = 1:num_bands
            selected = downsampling_sets{i}; 
            selected_signal=zeros(G.N,1);
            selected_signal(selected)=1;

            plot_param.climits = [0,1];
            figure;
            gsp_plot_signal(G, selected_signal, plot_param);
            title('Selected Vertices'); % add band i to title
        end
    end
end

% mcsfb_analysis
% Analysis function (takes G, filters, downsampling sets, returns transform
% coefficients; perform filtering and then downsample; return cell; total coefficients should be N)
% just one function (gsp_filter already does the switching)


param.shifted_ends = shifted_ends;
param.jackson = 1;
param.order=80;

[analysis_coeffs,filter_coeffs] = mcsfb_analysis(G, f, filter_bank, downsampling_sets, param);

if param.plot_filters
   approx_filters=cell(num_bands,1);
   xx=0:.001:G.lmax;
   filter_sum=zeros(size(xx));
   for i=1:num_bands
       approx_filters{i}=@(x) gsp_cheby_eval(x,filter_coeffs(:,i),[0,G.lmax]);
       filter_sum=filter_sum+approx_filters{i}(xx);
   end
   plot_param.show_sum=1;
   gsp_plot_filter(G,approx_filters,plot_param);
   hold on;
   plot(xx,filter_sum,'r:','LineWidth',3);   
   title('Approximate Filters');
end

% plot coeffs by channel
if param.plot_analysis_coeffs
    for i=1:num_bands
        figure;
        chan_coeffs=zeros(G.N,1);
        chan_coeffs(downsampling_sets{i})=analysis_coeffs{i};
        gsp_plot_signal(G,chan_coeffs,param);
        caxis([-2.5,2.5]);
        view(0,90)
        title('Analysis Coeffs by Channel');
    end
end


% mcsfb_sythesis =

[f_reconstruct, reconstruction_banded] = mcsfb_sythesis(G, num_bands, downsampling_sets, analysis_coeffs, shifted_ends, weights_banded);

% plot reconstruction for each channel
for i=1:num_bands
    figure;
    gsp_plot_signal(G,reconstruction_banded{i}, param);
    caxis([-2.5,2.5]);
    view(0,90)
    title('Reconstruction by Channel');
end

% plot reconstruction
figure;
max_val=max(abs(f));
plot_param.climits = [-max_val,max_val];
plot_param.vertex_size = 100;
gsp_plot_signal(G, f_reconstruct, plot_param);
caxis([-2.5,2.5]);
view(0,90)
title('Reconstruction');

% plot reconstruction error
error=abs(f-f_reconstruct);
plot_param.climits = [0, max(error)];
plot_param.vertex_size = 100;
figure;
gsp_plot_signal(G, error, plot_param);
caxis([-2.5,2.5]);
view(0,90)
title('Reconstruction Error');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dictionary atom
% exact_downsampling_partition

total_samples=length(vertcat(downsampling_sets{:}));
dict=zeros(G.N,total_samples);
for i=1:G.N
    del=gsp_delta(G,i);
    analysis_coeffs = mcsfb_analysis(G, del, filter_bank, downsampling_sets, param);
    dict(i,:)=vertcat(analysis_coeffs{:}); %sum(analysis_coeffs)';
end

atom_norms=sqrt(sum(dict.*dict,1))
G=gsp_compute_fourier_basis(G);
figure;
imagesc(G.U'*dict);
colorbar;
caxis([-.05,.05])

figure;
gsp_plot_signal_spectral(G,gsp_gft(G,dict(:,250)));

% normalized_dict=dict*diag(1./atom_norms);
% [~,ind]=sort(p,'ascend');
% sorted_dict=dict(:,ind);
% sorted_norm_dict=normalized_dict(:,ind);

% 
% % localization in graph Fourier domain? Yes
% figure;
% imagesc(G.U'*sorted_dict);
% colorbar;
% 
% 
% start=1;
% for i=1:length(subband_sizes)
%     figure;
%     plot(G.e,abs(G.U'*sorted_dict(:,start:start+subband_sizes(i)-1)));
%     hold on;
%     plot(G.e, mean(abs(G.U'*sorted_dict(:,start:start+subband_sizes(i)-1)),2),'k','LineWidth',3);
%     xlim([0,max(G.e)]);
%     set(gca,'FontSize',24);
%     start=start+subband_sizes(i);
%     xlabel('$\lambda$','Interpreter','latex','FontSize',24);
%     ylabel('$|\hat{\phi_i}(\lambda)|$','Interpreter','latex','FontSize',24);
% end
% 
% % localization in the vertex domain (normalized_atoms, sorted by frequency, with first row belonging to first two bands, and so forth)
% param.vertex_size=300;
% selected_inds=[1,160,325,627,1252];
% for i=1:length(selected_inds)
%     figure;
%     gsp_plot_signal(G,sorted_dict(:,selected_inds(i)),param);
%     %gsp_plot_signal(G,sorted_norm_dict(:,selected_inds(i)),param);
%     m=max(abs(sorted_dict(:,selected_inds(i))));
%     caxis([-m,m]);
%     view(0,90);
%     set(gca,'FontSize',24);
% end




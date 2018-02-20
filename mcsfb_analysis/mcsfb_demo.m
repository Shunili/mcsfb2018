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


% Params

% general
param.compute_full_eigen = 1;

% filter bank
num_bands = 4;
param.band_structure = 'method 2';
param.plot_filters = 1;

% downsampling
param.exact_downsampling_partition=0;
param.plot_downsampling_sets = 1;

% param.filter_type = 'approximate'; not used right now

if param.compute_full_eigen
    G = gsp_compute_fourier_basis(G);
else
    G = gsp_estimate_lmax(G);
end

% mcsfb_design_filter_bank
% Design filter bank (pass G, num_bands (M), parameters). Depend on whether
% you have all eigenvalues and what type of filter bank structure you want

[filter_bank,shifted_ends] = mcsfb_design_filter_bank(G,num_bands,param);

% plot filters
if param.plot_filters
    % add eigenvalues to plots?
    figure;
    plot_param.show_sum=0;
    gsp_plot_filter(G,filter_bank,plot_param);
end

% mcsfb_create_downsampling_sets
% Choose downsampling sets (either randomly or deterministically, depending
% on whether you have U; pass G and filter cell to this function). Output
% is a cell of length M where each entry has a vector with downsample
% vertices

downsampling_sets = mcsfb_create_downsampling_sets(G, filter_bank, shifted_ends, param);

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

analysis_coeffs = mcsfb_analysis(G, f, filter_bank, downsampling_sets, param);

% plot coeffs by channel
for i=1:num_bands
    figure;
    chan_coeffs=zeros(G.N,1);
    chan_coeffs(downsampling_sets{i})=analysis_coeffs{i};
    gsp_plot_signal(G,chan_coeffs,plot_param);
    caxis([-max(abs(chan_coeffs)),max(abs(chan_coeffs))]);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %check number of coeffs = G.N
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % Test if we found uniqueness sets
% Id=eye(G.N);
% parts=unique(subband_ids);
% num_part=length(parts);
% uniq_sets=ones(num_part,1);
% cond_nums=zeros(num_part,1);
% for i=1:num_part
%     T=[Id(:,p~=parts(i)),G.U(:,subband_ids==parts(i))];
%     if (size(T,1)~=G.N || size(T,2)~=G.N || (rank(T)~=G.N) )
%         uniq_sets(i)=0;
%     end
%     cond_nums(i)=cond(T);
% end
% cond_nums
% 
% if (sum(uniq_sets)==num_part)
%     display('Uniqueness partition successfully created');
% end
% 
% % test on one signal
% switch graph_type
%     case 'minnesota'
%         f=G.U(:,2)>0;
%     case 'sensor64'
%         f=G.U(:,1)+G.U(:,2)+G.U(:,5)+G.U(:,63);
%     case 'sensor500'
%         poly1a=-2*G.coords(:,1)+.5;
%         poly2a=G.coords(:,1).^2+G.coords(:,2).^2+.5;
%         f=zeros(G.N,1);
%         p1=(G.coords(:,2)>=(1-G.coords(:,1))) & (G.coords(:,2)<(1.5-G.coords(:,1)));
%         p2=(G.coords(:,2)<(0.6-G.coords(:,1)));
%         p3= p1 | p2;
%         f(~p3)=poly1a(~p3);
%         f(p3)=poly2a(p3);
%     otherwise
%         error('Unknown graph type');
% end
% coeffs=mcsfb(G,f,p,subband_ids);
% figure;
% gsp_plot_signal(G,f,param);
% figure;
% gsp_plot_signal(G,coeffs,param);
% 
% 
% % plot reconstruction for each channel
% cb=0;
% for i=1:num_part
%     figure;
%     vind=(p==parts(i));
%     sind=(subband_ids==parts(i));
%     if((sum(sind)>0) && (sum(vind)>0))
%         Q=G.U(vind,sind);
%         c=coeffs(vind);
%         alpha=Q\c;
%         interpolated_band_i = G.U(:,sind)*alpha;
%         cb=max(cb,max(abs(interpolated_band_i)));
%     else
%         interpolated_band_i=zeros(G.N,1);
%     end
%     gsp_plot_signal(G,interpolated_band_i,param);
%     caxis([-cb,cb]);
% end
% 
% 
% % examine resulting dictionary atoms (atoms in columns)
% dict=zeros(G.N);
% for i=1:G.N
%     del=gsp_delta(G,i);
%     c=mcsfb(G,del,p,subband_ids);
%     dict(i,:)=c';
% end
% 
% % full rank? Yes
% rank(dict)
% 
% % norms?
% atom_norms=sqrt(sum(dict.*dict,1))
% normalized_dict=dict*diag(1./atom_norms);
% [~,ind]=sort(p,'ascend');
% sorted_dict=dict(:,ind);
% sorted_norm_dict=normalized_dict(:,ind);
% 
% % orthogonal? Not quite
% % unnormalized correlation
% figure;
% imagesc(dict'*dict)
% colorbar;
% 
% % normalized correlation
% figure;
% Q=normalized_dict'*normalized_dict;
% imagesc(Q)
% caxis([-1,1]);
% colorbar;
% max(max(abs(Q-diag(diag(Q)))))
% 
% % sorted normalized correlation
% figure;
% Q=sorted_norm_dict'*sorted_norm_dict;
% imagesc(Q)
% caxis([-1,1]);
% colorbar;
% 
% % localization in graph Fourier domain? Yes
% figure;
% imagesc(G.U'*sorted_dict);
% colorbar;
% 
% 
% 
% % localization in the vertex domain
% 
% selected_inds=[1,32,65,126,252];
% for i=1:length(selected_inds)
%     figure;
%     param.vertex_size=50;
%     gsp_plot_signal(G,sorted_dict(:,selected_inds(i)),param);
% end
% 
% figure;
% param.vertex_size=50;
% for i=1:G.N
%     subplot(8,8,i);
%     gsp_plot_signal(G,dict(:,i),param);
%     caxis([-1,1]);
%     colorbar('off');
% end
% 
% % localization in the vertex domain (normalized_atoms, sorted by frequency, with first row belonging to first two bands, and so forth)
% figure;
% param.vertex_size=50;
% for i=1:G.N
%     subplot(8,8,i);
%     gsp_plot_signal(G,sorted_norm_dict(:,i),param);
%     caxis([-1,1]);
%     colorbar('off');
% end


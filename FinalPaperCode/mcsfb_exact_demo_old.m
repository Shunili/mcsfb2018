clear all;
close all;

rand('seed',0);
randn('seed',0);

% Initialize graph

% mcsfb_design_filter_bank
% Design filter bank (pass G, num_bands (M), parameters). Depend on whether
% you have all eigenvalues and what type of filter bank structure you want

% mcsfb_create_downsampling_sets
% Choose downsampling sets (either randomly or deterministically, depending
% on whether you have U; pass G and filter cell to this function). Output
% is a cell of length M where each entry has a vector with downsample
% vertices

% mcsfb_analysis
% Analysis function (takes G, filters, downsampling sets, returns transform
% coefficients; perform filtering and then downsample; return cell; total coefficients should be N)
% just one function (gsp_filter already does the switching)

% Just initialize graph here, not subbands
% Initialize graph and subbands 
graph_type='sensor64';
switch graph_type
    case 'minnesota'
        G=gsp_minnesota(1);
        subband_ids=[ones(165,1);2*ones(165,1);3*ones(330,1);4*ones(661,1);5*ones(1321,1)];
    case 'sensor64'
        G=gsp_sensor(64);
        subband_ids=[ones(4,1);2*ones(4,1);3*ones(8,1);4*ones(16,1);5*ones(32,1)];
    case 'sensor500'
        G=gsp_david_sensor_network(500);
        subband_ids=[ones(31,1);2*ones(31,1);3*ones(63,1);4*ones(125,1);5*ones(250,1)];
    otherwise
        error('Unknown graph type');
end

G=gsp_compute_fourier_basis(G);
plot_filters=0;

% Plot filters
if plot_filters
    num_filters=5;
    ideal_filter_bank=cell(num_filters,1);
    ideal_filter_bank{1}=@(x) ((0<=x) & (x<=G.e(31)));
    ideal_filter_bank{2}=@(x) ((G.e(32)<=x) & (x<=G.e(62)));
    ideal_filter_bank{3}=@(x) ((G.e(63)<=x) & (x<=G.e(125)));
    ideal_filter_bank{4}=@(x) ((G.e(126)<=x) & (x<=G.e(250)));
    ideal_filter_bank{5}=@(x) ((G.e(251)<=x) & (x<=G.e(500)));
    xx=0:.001:15;
    filter_data=zeros(length(xx),num_filters);
    for i=1:num_filters
        filter_data(:,i)=ideal_filter_bank{i}(xx);
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

%     ax.XTick = G.e;
%     set(gca,'XTickLabel','');
%     set(ax,'box','off','color','none')
%     % create new, empty axes with box but without ticks
%     bb = axes('Position',get(ax,'Position'),'box','on','xtick',[],'ytick',[]);
%     % set original axes as active
%     axes(ax)
%     % link axes in case of zooming
%     linkaxes([ax bb])
end

% Compute uniqueness set partition
p=uniqueness_set_partition_new(G,subband_ids);
figure;
param.vertex_size=100;
gsp_plot_signal(G,p,param);
hcb=colorbar;
colormap(prism(5));
caxis([.5,5.5]);
set(hcb,'YTick',1:5);
set(hcb,'FontSize',24);


% Test if we found uniqueness sets
Id=eye(G.N);
parts=unique(subband_ids);
num_part=length(parts);
uniq_sets=ones(num_part,1);
cond_nums=zeros(num_part,1);
for i=1:num_part
    T=[Id(:,p~=parts(i)),G.U(:,subband_ids==parts(i))];
    if (size(T,1)~=G.N || size(T,2)~=G.N || (rank(T)~=G.N) )
        uniq_sets(i)=0;
    end
    cond_nums(i)=cond(T);
end
cond_nums

if (sum(uniq_sets)==num_part)
    display('Uniqueness partition successfully created');
end

% test on one signal
switch graph_type
    case 'minnesota'
        f=G.U(:,2)>0;
    case 'sensor64'
        f=G.U(:,1)+G.U(:,2)+G.U(:,5)+G.U(:,63);
    case 'sensor500'
        poly1a=-2*G.coords(:,1)+.5;
        poly2a=G.coords(:,1).^2+G.coords(:,2).^2+.5;
        f=zeros(G.N,1);
        p1=(G.coords(:,2)>=(1-G.coords(:,1))) & (G.coords(:,2)<(1.5-G.coords(:,1)));
        p2=(G.coords(:,2)<(0.6-G.coords(:,1)));
        p3= p1 | p2;
        f(~p3)=poly1a(~p3);
        f(p3)=poly2a(p3);
    otherwise
        error('Unknown graph type');
end
coeffs=mcsfb(G,f,p,subband_ids);
figure;
gsp_plot_signal(G,f,param);
figure;
gsp_plot_signal(G,coeffs,param);

% plot coeffs by channel
for i=1:num_part
    figure;
    chan_coeffs=zeros(G.N,1);
    chan_coeffs(p==parts(i))=coeffs(p==parts(i));
    gsp_plot_signal(G,chan_coeffs,param);
    caxis([-max(abs(coeffs)),max(abs(coeffs))]);
end

% plot reconstruction for each channel
cb=0;
for i=1:num_part
    figure;
    vind=(p==parts(i));
    sind=(subband_ids==parts(i));
    if((sum(sind)>0) && (sum(vind)>0))
        Q=G.U(vind,sind);
        c=coeffs(vind);
        alpha=Q\c;
        interpolated_band_i = G.U(:,sind)*alpha;
        cb=max(cb,max(abs(interpolated_band_i)));
    else
        interpolated_band_i=zeros(G.N,1);
    end
    gsp_plot_signal(G,interpolated_band_i,param);
    caxis([-cb,cb]);
end


% examine resulting dictionary atoms (atoms in columns)
dict=zeros(G.N);
for i=1:G.N
    del=gsp_delta(G,i);
    c=mcsfb(G,del,p,subband_ids);
    dict(i,:)=c';
end

% full rank? Yes
rank(dict)

% norms?
atom_norms=sqrt(sum(dict.*dict,1))
normalized_dict=dict*diag(1./atom_norms);
[~,ind]=sort(p,'ascend');
sorted_dict=dict(:,ind);
sorted_norm_dict=normalized_dict(:,ind);

% orthogonal? Not quite
% unnormalized correlation
figure;
imagesc(dict'*dict)
colorbar;

% normalized correlation
figure;
Q=normalized_dict'*normalized_dict;
imagesc(Q)
caxis([-1,1]);
colorbar;
max(max(abs(Q-diag(diag(Q)))))

% sorted normalized correlation
figure;
Q=sorted_norm_dict'*sorted_norm_dict;
imagesc(Q)
caxis([-1,1]);
colorbar;

% localization in graph Fourier domain? Yes
figure;
imagesc(G.U'*sorted_dict);
colorbar;



% localization in the vertex domain

selected_inds=[1,32,65,126,252];
for i=1:length(selected_inds)
    figure;
    param.vertex_size=50;
    gsp_plot_signal(G,sorted_dict(:,selected_inds(i)),param);
end

figure;
param.vertex_size=50;
for i=1:G.N
    subplot(8,8,i);
    gsp_plot_signal(G,dict(:,i),param);
    caxis([-1,1]);
    colorbar('off');
end

% localization in the vertex domain (normalized_atoms, sorted by frequency, with first row belonging to first two bands, and so forth)
figure;
param.vertex_size=50;
for i=1:G.N
    subplot(8,8,i);
    gsp_plot_signal(G,sorted_norm_dict(:,i),param);
    caxis([-1,1]);
    colorbar('off');
end


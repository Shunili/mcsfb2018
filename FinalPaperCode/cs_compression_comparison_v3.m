% Code for numerical experiments in Section 5 fo Y. Jin and D. I Shuman, 
%  "An M-channel critically sampled filter bank for signals on graphs," 
%  ICASSP, 2017.

% This code requires the Graph Signal Processing Toolbox, which is
% available at https://lts2.epfl.ch/gsp/

% The spatial graph wavelet dictionary requires the MATLAB Boost Graph Library 
% (MatlabBGL) library from David Gleich, which is available at
% http://www.cs.purdue.edu/homes/dgleich/packages/matlab_bgl/index.html

% The QMF filterbank transform requires code from the USC Signal Transformation, 
% Analysis and Compression Group's demo package available 
% at http://biron.usc.edu/wiki/index.php/Graph_Filterbanks

% The diffusion wavelet transform requires Mauro Maggioni's diffusion 
% wavelet and diffusion geometry code, available at 
% http://www.math.duke.edu/~mauro/code.html

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures 11-13: Compression examples
close all;
clear all;
set(0,'DefaultFigureWindowStyle','docked');
rand('seed',40);
randn('seed',30);
 
% Can opt to not use these comparison transforms if the necessary toolboxes
% have not been installed. See above for links to the required codes
use_ckwt=0; % needs MatlabBGL 
use_qmf=1; % needs QMF filter bank code
use_diff_wav=1; % needs diffusion wavelet code
use_pyramid=0;
use_sgwt=0;

%% Create graph
example='pwSmoothSensor';

switch example
    case 'pwSmoothSensor'
        graph_param.N=500;
        G=gsp_david_sensor_network(graph_param.N);
            
    case 'smoothSensor'
        graph_param.N=500;
        G=gsp_david_sensor_network(graph_param.N);

    case 'localizedSensor'
        graph_param.N=500;
        G=gsp_david_sensor_network(graph_param.N);
        
    otherwise
        error('unknown example type');
end
G=gsp_compute_fourier_basis(G);

if use_pyramid
    %% Graph multiresolution
    num_levels=3;
    multiresolution_param.compute_full_eigen=1;
    multiresolution_param.sparsify=1;
    multiresolution_param.sparsify_epsilon=.3;
    Gs=gsp_graph_multiresolution(G,num_levels,multiresolution_param);

    makeExtraPlots=0;

    if makeExtraPlots
        for i=1:num_levels+1
            figure;
            gsp_plot_graph(Gs{i})
        end
    end
end

%% Create signal(s)
switch example
    case 'pwSmoothSensor'
        poly1a=-2*G.coords(:,1)+.5;
        poly2a=G.coords(:,1).^2+G.coords(:,2).^2+.5;
        signal=zeros(G.N,1);
        p1=(G.coords(:,2)>=(1-G.coords(:,1))) & (G.coords(:,2)<(1.5-G.coords(:,1)));
        p2=(G.coords(:,2)<(0.6-G.coords(:,1)));
        p= p1 | p2;
        signal(~p)=poly1a(~p);
        signal(p)=poly2a(p);
    
    case 'smoothSensor'
        signal=rand(G.N,1)-.5;
        tau=.5;
        filter=@(x) exp(-tau*x);
        signal=gsp_filter(G,filter,signal);
        
    case 'localizedSensor'
        % Generate ball radii
        radius=4;
        center_vertex=[13,444,22]; 
        num_balls=length(center_vertex);
        in_ball=zeros(G.N,num_balls);
        for i=1:num_balls
            in_ball(:,i)=sign((G.A^radius)*gsp_delta(G.N,center_vertex(i)));
        end

        % Generate filters
        num_filters=16;
        uniform_filters=gsp_design_half_cosine(G,num_filters);

        % Plot only the filters of interest
        filt_inds=[2,6,10,14];
        figure;
        gsp_plot_filter(G,uniform_filters(filt_inds));

        % Generate signal
        signal=zeros(G.N,1);
        actual_clusters=zeros(G.N,1);
        for i=1:num_balls+1
            signal_part=rand(G.N,1);
            signal_part=gsp_filter(G,uniform_filters{filt_inds(i)},signal_part);
            if i<=num_balls
                signal(in_ball(:,i)==1)=signal_part(in_ball(:,i)==1);
                actual_clusters(in_ball(:,i)==1)=i;
            else
                in_no_ball=max(0,1-sum(in_ball,2));
                signal(in_no_ball==1)=signal_part(in_no_ball==1);
                actual_clusters(in_no_ball==1)=i;
            end
        end
        
    case 'flickr'
        all_signals=signals(:,sum(signals)>0);
        signal=all_signals(:,26);
        num_signals=size(all_signals,2);
        
    otherwise
        error('unknown example type');
end

%% Plot the signal

% Figure (a): Plot the signal in the vertex domain
cmin=min(0,1.01*min(signal));
cmax=max(0,1.01*max(signal));
plot_param.climits=[cmin,cmax];
plot_param.vertex_size=100;
figure;
gsp_plot_signal(G,signal,plot_param);
%title('Signal','FontSize',24)
set(gca,'FontSize',24);

% Figure (b): Plot the signal in the graph spectral domain
figure;
signal_param.plot_abs=1;
gsp_plot_signal_spectral(G,gsp_gft(G,signal),signal_param);
xlabel('$\lambda$','Interpreter','latex','FontSize',24);
yy=ylabel('$|\hat{f}(\lambda)|~$','Interpreter','latex','FontSize',24);
%xlabel('$\lambda_{\ell}$','Interpreter','latex','FontSize',24);
%yy=ylabel('$|\hat{f}(\lambda_{\ell})|~$','Interpreter','latex','FontSize',24);
%set(gca,'FontSize',24)

%% Create MCSFB dictionary

% Design filter bank
num_bands = 5;
param.band_structure = 0;
param.plot_filters = 1;
[filter_bank,shifted_ends] = mcsfb_design_filter_bank(G,num_bands,param);

param.exact_downsampling_partition=1;
[downsampling_sets, weights_banded] = mcsfb_create_downsampling_sets(G, filter_bank, shifted_ends, param);

mcsfb_dict=zeros(G.N); % put atoms in rows
%allocated=0;
for i=1:num_bands
    %new_allocation=length(downsampling_sets{i});
    hmL=G.U*diag(filter_bank{i}(G.e))*G.U';
    mcsfb_dict(downsampling_sets{i},:)=hmL(downsampling_sets{i},:);
    %mcsfb_dict((allocated+1):(allocated+new_allocation),:)=hmL(downsampling_sets{i},:);
    %allocated=allocated+new_allocation;
end
% % have to change this for each example
% %subband_ids=[ones(31,1);2*ones(31,1);3*ones(63,1);4*ones(125,1);5*ones(250,1)];
% %subband_ids=[ones(7,1);2*ones(8,1);3*ones(16,1);4*ones(31,1);5*ones(63,1);6*ones(125,1);7*ones(250,1)];
% subband_ids=[ones(31,1);2*ones(31,1);3*ones(63,1);4*ones(125,1);5*ones(250,1)];
% vertex_partition=uniqueness_set_partition_new(G,subband_ids);
% %vertex_partition=subband_ids(randperm(G.N)); try uniform random sampling
% mcsfb_dict=zeros(G.N);
% for i=1:G.N
%     del=gsp_delta(G,i);
%     c=mcsfb(G,del,vertex_partition,subband_ids);
%     mcsfb_dict(i,:)=c';
% end
% % use extra step to normalize dict?
% %mcsfb_dict=mcsfb_dict*diag(1./sqrt(sum(mcsfb_dict.*mcsfb_dict,1)));
% mcsfb_dict=mcsfb_dict'; % rearranges dict so atoms are in rows

%% Create SGWT dictionary
% This is the spectral graph wavelet dictionary from D. K. Hammond, P.
% Vandergheynst, and R. Gribonval, "Wavelets on graphs via spectral graph
% theory," Appl. Comput. Harmon. Anal., vol. 30, no. 2, pp. 129-150, Mar.
% 2011.

if use_sgwt
    % Design filters for transform
    Nscales=6;
    sgwt_filters=gsp_design_abspline(G, Nscales);

    % Form SGWT dictionary
    sgwt_dict=zeros(Nscales*G.N,G.N);
    for i=1:Nscales
        sgwt_dict(G.N*(i-1)+1:G.N*i,:)=G.U*diag(sgwt_filters{i}(G.e))*G.U';
    end
end

%% Create QMF dictionary
% The QMF filterbank dictionary is based on the paper: S. K. Narang and 
% A.Ortega, "Perfect reconstruction two-channel wavelet filter banks for 
% graph structured data," IEEE Trans. Signal Process., vol. 60, no. 6, 
% pp. 2786-2799, Jun. 2012.

if use_qmf
    cheb_approx=40;
    qmf_dict=zeros(G.N,G.N); % atoms are the rows
    
    if strcmp(example,'flickr')
        qmf_param.color_method='dsatur'; % the BSC coloring method is too slow for this one
        [qmf_1,qmf_param]=usc_qmf_analysis(G,gsp_delta(G.N,1),cheb_approx,qmf_param);
    else
        [qmf_1,qmf_param]=usc_qmf_analysis(G,gsp_delta(G.N,1),cheb_approx);
    end
    
    qmf_dict(:,1)=qmf_1;
    for i=2:G.N
        qmf_dict(:,i)=usc_qmf_analysis(G,gsp_delta(G.N,i),cheb_approx,qmf_param);
    end
end

%% Create CKWT dictionary
% This is the spatial graph wavelet dictionary of M. Crovella and E.
% Kolaczyk, "Graph wavelets for spatial traffic analysis," in Proc. IEEE
% INFOCOM, vol. 3, Mar. 2003, pp. 1848-1857.
% Note: as designed the wavelets are all zero mean, so the analysis
% operator is not full rank. Therefore, we've added N Kronecker deltas,
% centered at each of the N vertices, to the dictionary.

if use_ckwt
    distances=zeros(G.N,G.N);
    for i=1:G.N
        distances(:,i)=shortest_paths(double(G.A),i); % shortest_paths function from MatlabBGL library
    end

    Nscales=5; % j runs from 1 to Nscales
    j_max=Nscales;
    rescale_factor=4; % how much to stretch the mexican_hat_wavelet
    wavelet_weights=zeros(Nscales,Nscales+1); %psi_jh for j running from 1 to j_max and h running from 0 to j
    for j=1:j_max
        for h=0:j
            wavelet_weights(j,h+1)=(j+1)*quad(@(x)mexican_hat_renormalized(x,rescale_factor),h/(j+1),(h+1)/(j+1));
        end
    end

    ckwt_dict=zeros((Nscales+1)*G.N,G.N); % atoms in the rows again
    for k=1:G.N % center vertex
        for j=1:j_max
            Cjk=0;
            for h=0:j
                Ikh=(distances(k,:)==h);
                ckwt_dict((j-1)*G.N+k,:)=ckwt_dict((j-1)*G.N+k,:)+Ikh*wavelet_weights(j,h+1)/sum(Ikh);
                Cjk=Cjk+wavelet_weights(j,h+1)^2/sum(Ikh);
            end
            Cjk=Cjk^(-.5);
            ckwt_dict((j-1)*G.N+k,:)=Cjk*ckwt_dict((j-1)*G.N+k,:);
        end
    end
    ckwt_dict(Nscales*G.N+1:end,:)=eye(G.N); % add all deltas to dictionary
end

%% Create diffusion wavelet dictionary
% This is the diffusion wavelet dictionary of R. R. Coifman and M.
% Maggioni, "Diffusion wavelets," Appl. Comput. Harmon. Anal., vol. 21, no.
% 1, pp. 53-94, Jul. 2006.

if use_diff_wav
    tau=1; % lower taus seem to work better for this graph. performance improves moving tau from 1 down to .5 to .2 to .1
    heat_kernel=@(x) exp(-tau*x); % change this -> j not explicitly set. should it be j_max?
    T=G.U*diag(heat_kernel(G.e))*G.U';
    j_max=6;
    Tree = DWPTree (T, j_max, 1e-4, ...
                    struct('Wavelets', true, 'OpThreshold', 1e-2, ...
                    'GSOptions', struct('StopDensity', 10, 'Threshold', 1e-3), ...
                    'Symm', true));
    % place atoms in the columns
    diff_wav_dict=zeros(G.N,G.N);
    num_dw_lev=size(Tree,1);
    num_scaling=size(Tree{num_dw_lev,1}.ExtBasis,2);
    diff_wav_dict(:,1:num_scaling)=Tree{num_dw_lev,1}.ExtBasis;
    start=num_scaling+1;
    for i=num_dw_lev:-1:1
        num_wav=size(Tree{i,2}.ExtBasis,2);
        diff_wav_dict(:,start:start+num_wav-1)=Tree{i,2}.ExtBasis;
        start=start+num_wav;
    end
    diff_wav_dict=diff_wav_dict'; % change atoms to rows
end

%% Create pyramid dictionary
if use_pyramid
    analysis_param.h_filters=cell(num_levels,1);
    for i=1:num_levels
        meyer_param.t=4/(3*Gs{i}.lmax);
        meyer_filt=gsp_design_meyer(Gs{i},2,meyer_param);
        analysis_param.h_filters{i}=meyer_filt{1};
    end

    analysis_param.regularize_epsilon=.005; 
    syntehsis_param.h_filters=analysis_param.h_filters;
    synthesis_param.regularize_epsilon=analysis_param.regularize_epsilon;

    % Form pyramid dictionary
    num_pyr_atoms=0;
    for i=1:num_levels+1
        num_pyr_atoms=num_pyr_atoms+Gs{i}.N;
    end
    pyramid_dict=zeros(num_pyr_atoms,G.N);
    for i=1:G.N
        [coarse_approximations,prediction_errors]=gsp_pyramid_analysis(Gs,gsp_delta(G.N,i),num_levels,analysis_param);
        pyramid_dict(:,i)= [coarse_approximations{num_levels+1} ; vertcat(prediction_errors{:})];
    end
end

%% Compute all transform coefficients and produce rate/distortion curves

% Compute all transform coefficients
mcsfb_coeffs=mcsfb_dict*signal;
if use_sgwt
    sgwt_coeffs=sgwt_dict*signal;
end
gft_coeffs=gsp_gft(G,signal);
if use_pyramid
    pyramid_coeffs=pyramid_dict*signal;
end
if use_ckwt
    ckwt_coeffs=ckwt_dict*signal;
end
if use_qmf
    qmf_coeffs=qmf_dict*signal;
end
if use_diff_wav
    diff_wav_coeffs=diff_wav_dict*signal;
end

% Figure (c): Plot sorted coefficients and compare to pyramid transform
mcsfb_coeffs_sorted=sort(abs(mcsfb_coeffs),'descend');
if use_sgwt
    sgwt_coeffs_sorted=sort(abs(sgwt_coeffs),'descend');
end
if use_pyramid
    pyramid_coeffs_sorted=sort(abs(pyramid_coeffs),'descend');
end
gft_coeffs_sorted=sort(abs(gft_coeffs),'descend');
deltas_coeffs_sorted=sort(abs(signal),'descend');
if use_ckwt
    ckwt_coeffs_sorted=sort(abs(ckwt_coeffs),'descend');
end
if use_qmf
    qmf_coeffs_sorted=sort(abs(qmf_coeffs),'descend');
end
if use_diff_wav
    diff_wav_coeffs_sorted=sort(abs(diff_wav_coeffs),'descend');
end

figure;
hold on;
plot(1:size(mcsfb_dict,1),mcsfb_coeffs_sorted/mcsfb_coeffs_sorted(1),'-r','LineWidth',3);

if use_pyramid
    plot(1:size(pyramid_dict,1),pyramid_coeffs_sorted/pyramid_coeffs_sorted(1),'-r','LineWidth',3);
end
plot(1:G.N,gft_coeffs_sorted/gft_coeffs_sorted(1),'-g','LineWidth',3);
plot(1:G.N,deltas_coeffs_sorted/deltas_coeffs_sorted(1),'-','color',[.35,.35,.35], 'LineWidth',3);
if use_sgwt
    plot(1:size(sgwt_dict,1),sgwt_coeffs_sorted/sgwt_coeffs_sorted(1),':k','LineWidth',3);
end
if use_ckwt
    plot(1:size(ckwt_dict,1),ckwt_coeffs_sorted/ckwt_coeffs_sorted(1),'-.m','LineWidth',3);
end
if use_qmf
    plot(1:G.N,qmf_coeffs_sorted/qmf_coeffs_sorted(1),'-.b','LineWidth',3);
end
if use_diff_wav
    plot(1:G.N,diff_wav_coeffs_sorted/diff_wav_coeffs_sorted(1),':c','LineWidth',3);
end
%legend('Pyramid','GFT','Deltas','SGWT','CKWT','QMF','DiffWav');
legend('MCSFB','GFT','Deltas','QMF','DiffWav');
set(gca,'FontSize',24)
box on;
ylabel('Normalized Magnitude');
xlabel('Sorted Coefficient Index');

% Figure (d): Produce rate/distortion curves
if use_sgwt
%        sgwt_errors=compute_all_hard_thresholds(signal,sgwt_dict,sgwt_coeffs);  
    sgwt_errors=compute_all_omp(signal,sgwt_dict);

end
%mcsfb_errors=mcsfb_compute_all_hard_thresholds(G,mcsfb_coeffs,vertex_partition,subband_ids,signal);
%mcsfb_errors=compute_all_hard_thresholds(signal,mcsfb_dict,mcsfb_coeffs);
mcsfb_errors=compute_all_omp(signal,mcsfb_dict);
gft_errors=compute_all_omp(signal,G.U',gft_coeffs);
delta_errors=compute_all_omp(signal,eye(G.N),signal);
if use_pyramid
    %pyramid_errors=compute_all_hard_thresholds(signal,pyramid_dict,pyramid_coeffs);
    pyramid_errors=compute_all_omp(signal,pyramid_dict);
end
if use_ckwt
    ckwt_errors=compute_all_omp(signal,ckwt_dict,ckwt_coeffs);
end
if use_qmf
    qmf_errors=compute_all_omp(signal,qmf_dict,qmf_coeffs);
end
if use_diff_wav
    diff_wav_errors=compute_all_omp(signal,diff_wav_dict,diff_wav_coeffs);
end

figure;
hold on;
plot(1:G.N,mcsfb_errors.^2,'-r','LineWidth',3);
plot(1:G.N,gft_errors.^2,'-g','LineWidth',3);
plot(1:G.N,delta_errors.^2,'-','color',[.35,.35,.35], 'LineWidth',3);
if use_pyramid
    plot(1:G.N,pyramid_errors(1:G.N).^2,'-y','LineWidth',3);
end
if use_sgwt
    plot(1:G.N,sgwt_errors(1:G.N).^2,':k','LineWidth',3);
end
if use_ckwt
    plot(1:G.N,ckwt_errors.^2,'-.m','LineWidth',3);
end
if use_qmf
    plot(1:G.N,qmf_errors.^2,'-.b','LineWidth',3);
end
if use_diff_wav
    plot(1:G.N,diff_wav_errors.^2,':c','LineWidth',3);
end
legend('MCSFB','GFT','Deltas','QMF','DiffWav');
%legend('MCSFB','GFT','Deltas','Pyramid','SGWT','QMF','DiffWav');
%legend('Pyramid','GFT','Deltas','SGWT','CKWT','QMF','DiffWav');
set(gca,'FontSize',24)
box on;
ylabel('NMSE');
xlabel('Number of Coefficients Used in Reconstruction');
axis([0 500 0 1]);

% plot on log axis
cap=400;

figure;
semilogy(1:cap,mcsfb_errors(1:cap).^2,'-r','LineWidth',3);
hold on;
plot(1:cap,gft_errors(1:cap).^2,'-g','LineWidth',3);
plot(1:cap,delta_errors(1:cap).^2,'-','color',[.35,.35,.35], 'LineWidth',3);
if use_pyramid
    plot(1:cap,pyramid_errors(1:cap).^2,'-y','LineWidth',3);
end
if use_sgwt
    plot(1:cap,sgwt_errors(1:cap).^2,':k','LineWidth',3);
end
if use_ckwt
    plot(1:cap,ckwt_errors(1:cap).^2,'-.m','LineWidth',3);
end
if use_qmf
    plot(1:cap,qmf_errors(1:cap).^2,'-.b','LineWidth',3);
end
if use_diff_wav
    plot(1:cap,diff_wav_errors(1:cap).^2,':c','LineWidth',3);
end
legend('MCSFB','GFT','Deltas','QMF','DiffWav','Location','SouthWest');
set(gca,'FontSize',24)
box on;
ylabel('NMSE');
xlabel('Number of Coefficients Used in Reconstruction');
grid on;

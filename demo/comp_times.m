close all;
clear all;

rand('seed',0);
randn('seed',0);

% Main parameters to explore
adapted=1; % downsampling sets and number of measurements adapted to signal
exact = 0;
scenarioA = 0;
method='gft'; % 'mcsfb' or 'diffusion' or 'qmf' of 'gft'
graph='net25';

% Other parameters
num_bands = 5;
param.search_right_only=0;
param.replacement=0; % use 1 for large graphs for now as the without replacement method is too slow
synth_param.reg_filter=3; % reconstruction method: splines
synth_param.gamma=1; % surprisingly insensitive to this parameter % larger gamma puts more weight on matching samples; smaller gamma puts more weight on matching spectral content; convergence must faster for larger gamma, which makes sense since we are initializing it to the guesses
param.subtract_mean=1;
param.num_vec=30; % default of 30 for most reconstruction methods seems fine. Plays a more important role if we are trying to reconstruct from the approximated subspaces
make_plots=0;

% Graph and signal

switch graph
    case 'gnp'
        N=100;
        p=.3;
        G=gsp_erdos_renyi(N,p);
    case 'sensor'
        N=500;
        G=gsp_david_sensor_network(N);
        sig=randn(G.N,1);
        signal=G.L*sig;
        vs=80;
        plim=[-25,25];
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
        comm_param.Nc=100;
        N=25000;
        G=gsp_community(N,comm_param);
        signal=randn(G.N,1);
    case 'net25' 
%        load('/Users/davidshuman/Dropbox/Current Research Work/MCSFB/Shuni_Thesis/GitHub/mcsfb2018/net25_data/net25.mat');
        A=Problem.A;
        A = A - diag(diag(A)); 
        A(4228,6327) = 1;
        A(6327,4228) = 1;
        G=gsp_graph(A);
        signal=randn(G.N,1);
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
        [V,D]=eigs(G1.L+1e-11*speye(G1.N),1,'sm'); % make sure this works for other machines. compute total number in each and choose the large one
        large_component=(abs(V(:,1))>.000000001);
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

if exact
    adapted=0;
    param.exact_downsampling_partition=1;
elseif scenarioA 
    param.order=25; % used for density estimation and analysis filtering
    synth_param.order=param.order; 
    synth_param.pcgtol=1e-8; 
    synth_param.pcgmaxits=100;
else
    param.order=50; %100; % used for density estimation and analysis filtering
    synth_param.order=param.order; 
    synth_param.pcgtol=1e-10;  
    synth_param.pcgmaxits=250; 
end

initial_signal=signal;        
switch method
    case 'gft'
        tic
        G=gsp_compute_fourier_basis(G);
        fhat=gsp_gft(G,signal);
        setup_analysis_time=toc
        tic
        f_reconstruct=gsp_igft(G,fhat);
        synth_time=toc
    case 'qmf'
    cheb_approx=50;
    tic
    %if G.N>1000
     
    qmf_param.color_method='dsatur'; % the BSC coloring method is too slow for this one
    [qmf_1,qmf_param]=usc_qmf_analysis(G,signal,cheb_approx,qmf_param);
    %else
    %    [qmf_1,qmf_param]=qmf_analysis(G,signal,cheb_approx);
    %end
    setup_analysis_time=toc
    
    case 'diffusion'
        tic
        G=gsp_compute_fourier_basis(G);
        % diffusion operator
        tau=1;
        heat_kernel=@(x) exp(-tau*x);
        T=G.U*diag(heat_kernel(G.e))*G.U';
        Tree = DWPTree (T, num_bands-1, 1e-4, ...
                struct('Wavelets', true, 'OpThreshold', 1e-2, ...
                'GSOptions', struct('StopDensity', 1, 'Threshold', 1e-3), ...
                'Symm', true));
        setup_time=toc
        tic
        CoeffTree = DWCoeffs(Tree, signal);
        analysis_time=toc
        setup_analysis_time=setup_time+analysis_time
        tic
        f_reconstruct = DWRecon(Tree, DWWavelet(CoeffTree));
        synth_time=toc
        
    case 'mcsfb'
        if make_plots
            figure;
            param.vertex_size=vs;
            param.climits=plim;
            gsp_plot_signal(G,signal,param);
            set(gca,'FontSize',24);
            if strcmp(graph,'bunny')
                view(0,90);
            end
            title('Initial signal');
        end

        if param.subtract_mean
            signal_average=mean(signal)
            signal=signal-signal_average;
            if make_plots
                figure;
                param.vertex_size=vs;
                newlim=max(abs(signal));
                param.climits=[-newlim,newlim];
                gsp_plot_signal(G,signal,param);
                set(gca,'FontSize',24);
                if strcmp(graph,'bunny')
                    view(0,90);
                end
                title('Initial mean normalized signal');
            end
        end

        % design filter bank
        tic
        if exact
            G=gsp_compute_fourier_basis(G);
        else
            G = gsp_estimate_lmax(G);
        end
        param.band_structure = 0;
        param.spectrum_adapted=1;


        [filter_bank,shifted_ends,band_ends,G] = mcsfb_design_filter_bank(G,num_bands,param);
        filter_design_time=toc

        % create downsampling sets (and ensure critical sampling)
        tic
        param.jackson = 1;
        param.shifted_ends = shifted_ends;
        if adapted
            [param.signal_projections,filter_coeffs]=mcsfb_apply_filters(G,signal,filter_bank,param);
            param.adapt_weights=1;
            param.adapt_num_meas=1;  
        end
        param.target_samples=G.N-param.subtract_mean;
        [downsampling_sets, weights_banded] = mcsfb_create_downsampling_sets(G, filter_bank, shifted_ends, param);
        downsampling_selection_time=toc

        total_setup_time=filter_design_time+downsampling_selection_time

        % perform analysis
        tic
        if adapted % already computed filter coeffs
            [analysis_coeffs] = mcsfb_analysis(G, signal, filter_bank, downsampling_sets, param);
        else
            [analysis_coeffs,filter_coeffs] = mcsfb_analysis(G, signal, filter_bank, downsampling_sets, param);
        end
        analysis_time=toc

        total_setup_analysis_time=total_setup_time+analysis_time

        % synthesis
        tic
        [f_reconstruct, reconstruction_banded] = mcsfb_synthesis(G, num_bands, downsampling_sets, analysis_coeffs, shifted_ends, weights_banded, synth_param);
        if param.subtract_mean
            f_reconstruct=f_reconstruct+signal_average;
        end
        synth_time=toc

        % tic
        % f_reconstruct2 = mcsfb_synth2(G, downsampling_sets, filter_coeffs, analysis_coeffs, synth_param);
        % if param.subtract_mean
        %     f_reconstruct2=f_reconstruct2+signal_average;
        % end
        % synth2_time=toc
end

 % compute reconstruction error
error=f_reconstruct-initial_signal;
mse=sum(error.^2)/G.N
nmse=mse/sum(signal.^2)

% plot overall reconstruction and error
if make_plots
    param.climits=plim;
    figure;
    gsp_plot_signal(G,f_reconstruct,param);
    set(gca,'FontSize',24);
    title('Reconstruction');
    if strcmp(graph,'bunny')
        view(0,90);
    end

    error_param.vertex_size=vs;
    figure;
    gsp_plot_signal(G,abs(error),error_param);
    set(gca,'FontSize',24);
    title('Reconstruction error');
    if strcmp(graph,'bunny')
        view(0,90);
    end
    colormap(flipud(hot));
end
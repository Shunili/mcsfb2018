function downsampling_sets = mcsfb_create_downsampling_sets(G, filter_bank, shifted_ends, param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 
    num_bands=size(filter_bank,1);
    downsampling_sets=cell(num_bands,1);  
    
    if ~isfield(param,'exact_downsampling_partition')
        param.exact_downsampling_partition=0;
    end
    
    if param.exact_downsampling_partition
        
        if ~gsp_check_fourier(G)
            G=gsp_compute_fourier_basis(G); 
        end
        subband_ids = zeros(G.N, 1);
        for i=1:num_bands
            subband_ids(filter_bank{i}(G.e)~=0)=i; % if some eigenvalues are included in multiple filters, this code will put them in the last one 
        end
%         subband_ids(1:floor(G.N*(1/2)^(num_bands-1)),:) = ones(floor(G.N*(1/2)^(num_bands-1)), 1);
%          t = 2;
%         for i=num_bands-1:-1:0
%             subband_ids((floor(G.N*(1/2)^(i+1))+1):floor(G.N*(1/2)^i),:) = t*ones(floor(G.N*(1/2)^i)-(floor(G.N*(1/2)^(i+1))+1)+1,1);
%             t = t + 1;
%         end
        partition_ids = part_mat(G.U,subband_ids,param);
        for i=1:num_bands
            downsampling_sets{i}=find(partition_ids==i);
        end
        
    else
        
        exact=gsp_check_fourier(G);
        if ~exact
            G=gsp_spectrum_cdf_approx(G);
        end
        
        %update name
        num_its = ceil(2*log(G.N));
        
        for i = 1:num_bands
            
            h = filter_bank{i};

           % r = rand(G.N,1);
           % y = gsp_filter(G,h,r); 
            
            up_limit = shifted_ends(i+1);
            low_limit = shifted_ends(i);
            
            % find approximate number of eigenvalues in each band
            if exact
                %extra_samps=0
                nb_meas = sum(h(G.e)); %+extra_samps; % m: num eigenvalues in band, will need to estimate if don't have exact eigenvalues
            else
                nb_meas = floor((G.spectrum_cdf_approx(up_limit)-G.spectrum_cdf_approx(low_limit))*G.N);
            end

            [weights, ~] = compute_sampling_weights(G,num_its,h);
            [~, selected] = build_sampling_matrix(G, weights, nb_meas);

            downsampling_sets{i} = selected;
        end
    end
 end
  


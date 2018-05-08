function [downsampling_sets, second_output] = mcsfb_create_downsampling_sets(G, filter_bank, shifted_ends, param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 
    num_bands=size(filter_bank,1);
    downsampling_sets=cell(num_bands,1); 
    second_output = cell(num_bands,1); %weights banded if not exact, eigenvector ids for each band for exact
    
    if ~isfield(param,'exact_downsampling_partition')
        param.exact_downsampling_partition=0;
    end
    
    if ~isfield(param,'prenormalize_num_meas')
        param.prenormalize_num_meas=1;
    end
    
%     if ~isfield(param,'alpha')
%         param.alpha=.25;
%     end
    
    if ~isfield(param,'fixed_total_meas')
        param.fixed_total_meas=1;
    end
    
    if ~isfield(param,'subtract_mean')
        param.subtract_mean=0;
    end
    
    if ~isfield(param,'target_samples')
        param.target_samples=G.N;
    end
    
    if ~isfield(param,'extra_low_factor')
        param.extra_low_factor=1;
    end
    
    if ~isfield(param,'adapt_weights')
        param.adapt_weights=0;
    else
        if param.adapt_weights
            if ~isfield(param,'signal_projections')
                error('param.signal_projections must be defined to use signal-adapted weights');
            end
        end
    end
    
    if ~isfield(param,'adapt_num_meas')
        param.adapt_num_meas=0;
    else
        if param.adapt_num_meas
            if ~isfield(param,'signal_projections')
                error('param.signal_projections must be defined to use signal-adapted number of measurements');
            end
        end
    end
    
    if param.exact_downsampling_partition
        
        if ~gsp_check_fourier(G)
            G=gsp_compute_fourier_basis(G); 
        end
        subband_ids = zeros(G.N, 1);
        for i=1:num_bands
            subband_ids(filter_bank{i}(G.e)~=0)=i;% if some eigenvalues are included in multiple filters, this code will put them in the last one 
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
            second_output{i}=(filter_bank{i}(G.e)~=0);
        end
        
    else
        weights = cell(num_bands,1);
        norm_Uk = cell(num_bands,1);
        nb_meas = zeros(num_bands,1);
        
        exact=gsp_check_fourier(G);
        

        for i = 1:num_bands % compute both the weights for each band and an initial number of measurements for each band
            
            if exact
                %extra_samps=0
                h = filter_bank{i};
                nb_meas(i) = sum(h(G.e)); %+extra_samps; % m: num eigenvalues in band, will need to estimate if don't have exact eigenvalues
                norm_Uk{i}= sum(G.U.^2, 2);
                weights{i}=norm_Uk/sum(norm_Uk);

            else
                %nb_meas = round((G.spectrum_cdf_approx(up_limit)-G.spectrum_cdf_approx(low_limit))*G.N);
                
                up_limit = shifted_ends(i+1);
                low_limit = shifted_ends(i);
                [~, jch] = gsp_jackson_cheby_coeff(low_limit,up_limit,[0,G.lmax], param.order);
                
                % TODO: add option to make the
                % number of random vectors higher than the number used to
                % compute the density
                r=gsp_cheby_opX(G,jch);
                nb_meas(i)=round(gsp_hutch(G,r));

                % [weights, ~] = compute_sampling_weights(G,num_its,h);
                norm_Uk{i}= sum(r.^2, 2);
                weights{i}=norm_Uk{i}/sum(norm_Uk{i});
                if param.adapt_weights && (param.subtract_mean || (i>1)) 
                    weights{i}=norm_Uk{i}.*abs(param.signal_projections(:,i));
                    weights{i}=weights{i}/sum(weights{i});
                end
            end
           
        end
        
        second_output=weights;
                
        if ~exact % adjust the number of measurements     

            if param.adapt_num_meas % adapted to signal energy (forced to have target number of samples)
                proj_norms=sqrt(sum(param.signal_projections.^2));
                %nb_meas=nb_meas.*(proj_norms'.^param.alpha);
                nb_meas=nb_meas.*(log(proj_norms'));
                param.prenormalize_num_meas=1;
                param.fixed_total_meas=1;
            end

            total_samples=sum(nb_meas); % count number initially allocated
            if param.prenormalize_num_meas
                nb_meas=round(param.target_samples/total_samples*nb_meas);
            end
            
            total_samples=sum(nb_meas); 
            if param.fixed_total_meas
                % give to the low, take from the high
                if total_samples>param.target_samples % eliminate from last band
                    extra=total_samples-param.target_samples;
                    nb_meas(num_bands)=nb_meas(num_bands)-extra;
                elseif total_samples<param.target_samples % add more to first band
                    nb_meas(1)=nb_meas(1)+param.target_samples-total_samples;
                end
            end
            

        end
  
        for i = 1:num_bands % perform the sampling
            [~, selected] = build_sampling_matrix(G, weights{i}, nb_meas(i), param);
            downsampling_sets{i} = selected;
        end
    end
 end
  


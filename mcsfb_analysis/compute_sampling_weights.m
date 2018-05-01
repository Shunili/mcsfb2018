function [weights, P_min_half] = compute_sampling_weights(G,L,h)

    %sqnorm_UR = zeros(G.N,L);
    param.grid_order=1000;
    param.order=30;
    %tally=zeros(G.N,1);
    
    Sig = randn(G.N,L)*1/sqrt(L);
    X = gsp_filter(G,h,Sig,param);
    norm_Uk= sum(X.^2, 2);
    weights=norm_Uk/sum(norm_Uk);
    
%     f=randn(G.N,L);
%     f2=gsp_filter(G,h,f,param);
%     
    
%     for num = 1:L
%         f=randn(G.N,1);
%         f2 = gsp_filter(G,h,f,param);
%         tally=tally+f2.^2;
%         % Random signal
% %        f = randn(G.N);
%         % filter
%  %       f2 = gsp_filter(G,h,f,param);
%         % store 
%  %       sqnorm_UR(:, num) = sum(f2.^2, 2);
%     end
%         weights=tally/sum(tally);
%         %cum_coh = mean(sqnorm_UR, 2);
%         %weights = cum_coh/sum(cum_coh(:));
        

        P_min_half = sparse(1:G.N, 1:G.N, 1./sqrt(weights), G.N, G.N);
end


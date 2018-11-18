function [f_reconstruct] = mcsfb_synth2(G, downsampling_sets, filter_coeffs, analysis_coeffs, param)

if nargin<5
   param = struct;
end

if ~isfield(param,'pcgtol')
    param.pcgtol=1e-14;
end

if ~isfield(param,'pcgmaxits')
    param.pcgmaxits=5000;
end

num_bands=size(filter_coeffs,2);
right_side=zeros(G.N,1);
for m=1:num_bands
    analysis_upsampled=zeros(G.N,1);
    analysis_upsampled(downsampling_sets{m})=analysis_coeffs{m};
    right_side=right_side+gsp_cheby_op(G,filter_coeffs(:,m),analysis_upsampled,param);
end

LHS=@(z) 0;
for m=1:num_bands
    MtM=sparse(downsampling_sets{m},downsampling_sets{m},ones(length(downsampling_sets{m}),1),G.N,G.N);
    LHS=@(z) LHS(z)+gsp_cheby_op(G,filter_coeffs(:,m),MtM*gsp_cheby_op(G,filter_coeffs(:,m),z,param),param);
end

% estimate diagonal of LHS (Bekas, 2005)
numer=zeros(G.N,1); 
den=zeros(G.N,1); 
if isfield(G,'X')
    num_vec=size(G.X,2);
else
    if ~isfield(param,'num_vec')
        num_vec = 30; 
    else
        num_vec = param.num_vec;
    end
    G.X=randn(G.N, num_vec);
end
for i=1:num_vec 
    numer=numer+G.X(:,i).*LHS(G.X(:,i)); 
    den=den+G.X(:,i).^2; 
end

diag_est=numer./den;
preconditioner=@(z) z./diag_est;
f_reconstruct=pcg(LHS,right_side,param.pcgtol,param.pcgmaxits,preconditioner);



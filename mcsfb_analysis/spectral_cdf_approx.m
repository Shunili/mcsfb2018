function [ cdf_approx, vals ] = spectral_cdf_approx( G , param)

if ~isfield(param, 'order')
    param.order = 30;
end

if ~isfield(param, 'num_vec')
    param.num_vec = 30;
end

if ~isfield(param, 'pts')
    param.num_pts = 25;
    param.pts=linspace(0,G.lmax,param.num_pts);
end

if ~isfield(param, 'num_pts')
    param.num_pts = length(param.pts);
end

if ~isfield(G,'lmax')
    G=gsp_estimate_lmax(G);
end

% replace num_pts by actual pts, currently using linear spacing
% put in one more parameter, actual pts
%pts=linspace(0,G.lmax,param.num_pts);
vals=zeros(param.num_pts,1);
vals(1)=1;
%TODO: force last one to be N?

jch=zeros(param.order+1,param.num_pts-1);
Sig = randn(G.N,param.num_vec);
    
for j=1:param.num_pts-1
    [~, jch(:,j)] = gsp_jackson_cheby_coeff(0,param.pts(j+1),[0,G.lmax], param.order);
end

%this can be more efficient by computing T_bar(L)f once and store them.
%update the chby_op to retuen T_bar(L)f
X = gsp_cheby_op(G, jch, Sig); %jch, each column is the coeffs for one filter

St=Sig';
for j=1:param.num_pts-1
    for i=1:param.num_vec
        vals(j+1)=vals(j+1)+St(i,:)*X((j-1)*G.N+1:j*G.N,i)/param.num_vec;
    end
end
vals=vals/G.N;
if(vals(1)>vals(2))
    vals(1)=vals(2);
end
cdf_approx = @(s) gsp_mono_cubic_warp_fn(param.pts',vals,s);

end


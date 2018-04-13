function [ z ] = mcsfb_reconstruct_band2( G , selected, values, lower, upper, weights, param)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin<7
   param = struct;
end

if ~isfield(param,'exact_downsampling_partition')
    param.exact_downsampling_partition=0;
end

num_selected=length(selected);

if ~param.exact_downsampling_partition
    if length(values) ~= num_selected
        error('Values should be defined only on selected vertices');
    end
    if length(weights) ~= num_selected
        error('Weights should be defined only on selected vertices');
    end
else
    weights=ones(num_selected,1);
end
    
if ~isfield(param,'reg_eps')
    %reg_eps=1/G.N;
    reg_eps=1/1;
else
    reg_eps=param.reg_eps;
end

if ~isfield(param,'pcgtol')
    param.pcgtol=1e-14;
end

if ~isfield(param,'pcgmaxits')
    param.pcgmaxits=5000;
end

if ~isfield(param,'order')
    order=20;
else
    order=param.order;
end

h = @(x) (x>=lower & x<upper);
reg_filter =@(x) 1./(h(x)+reg_eps)-1/(1+reg_eps);
wd = zeros(G.N,1);
wd(selected)=1./weights;
B=diag(wd);
right_side = zeros(G.N,1);
right_side(selected)=values./weights;


if (isfield(G,'U') && isfield(G,'e'))
    A=B+G.U*diag(reg_filter(G.e))*G.U'; 
    z=A\right_side;
else
    range=[0,G.lmax];
    grid_order=1000;
    [~, JCH]=gsp_jackson_cheby_coeff(lower, upper, range, order);
    reg_short=@(x) 1./(gsp_cheby_eval(x,JCH,range)+reg_eps);
    short_coeff=gsp_cheby_coeff(G,reg_short,order,grid_order);  % may need to use damping here
    %gte=@(x) gsp_cheby_op(G, JCH, x)+reg_eps*x;
    %LHS=@(z) B*z-(1/(1+reg_eps))*z+pcg(gte,z,1e-10,100);
    LHS=@(z) B*z-(1/(1+reg_eps))*z+gsp_cheby_op(G,short_coeff,z);
    z=pcg(LHS,right_side,1e-10,100);
end




end


function [ z ] = mcsfb_reconstruct_band2( G , selected, values, lower, upper, weights, param)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin<7
   param = struct;
end

num_selected=length(selected);
if length(values) ~= num_selected
    error('Values should be defined only on selected vertices');
end
if length(weights) ~= num_selected
    error('Weights should be defined only on selected vertices');
end

if ~isfield(param,'reg_eps')
    reg_eps=1/G.N;
    reg_eps=1/100;
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
    order=80;
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
    [~, JCH]=gsp_jackson_cheby_coeff(lower, upper, range, order);
    gte=@(x) gsp_cheby_op(G, JCH, x)+reg_eps*x;
    LHS=@(z) B*z-(1/(1+reg_eps))*z+pcg(gte,z);
    z=pcg(LHS,right_side,1e-8,100);
end




end


function [ z ] = mcsfb_reconstruct_band( G , selected, values, lower, upper, weights, param)
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

if ~isfield(param,'gamma')
    gamma=G.N;
else
    gamma=param.gamma;
end

if ~isfield(param,'pcgtol')
    param.pcgtol=1e-14;
end

if ~isfield(param,'pcgmaxits')
    param.pcgmaxits=5000;
end

if ~isfield(param,'order')
    order=50;
else
    order=param.order;
end

h = @(x) (x>=lower & x<upper);
reg_filter =@(x) 1-h(x);
wd = zeros(G.N,1);
wd(selected)=1./weights;
B=diag(wd);
right_side = zeros(G.N,1);
right_side(selected)=values./weights;

if (isfield(G,'U') && isfield(G,'e'))
    A=B+gamma*G.U*diag(reg_filter(G.e))*G.U'; 
    z=A\right_side;
else % use conjugate gradient to solve linear system of equations, need to be careful because Chebyshev polynomial approx can result in negative values, making A indefinite. Jackson coeffs should keep A positive definite
    if ~isfield(param,'precondition')
        precondition=1;
    else
        precondition=param.precondition;
    end
    range=[0,G.lmax];
    [~, JCH]=gsp_jackson_cheby_coeff(lower, upper, range, order);
    JCH_reg=zeros(size(JCH));
    JCH_reg(1)=2;
    JCH_reg=JCH_reg-JCH;
    reg_filter_eval= @(x)gsp_cheby_op(G,JCH_reg,x);
    afun = @(x) B*x+gamma*reg_filter_eval(x);
    if precondition
        preconditioner=@(z) z./(wd+gamma);
        initial_guess=zeros(G.N,1);
        initial_guess(selected)=values;
        z=pcg(afun,right_side,param.pcgtol,param.pcgmaxits,preconditioner,[],initial_guess);
    else
        z=pcg(afun,right_side,param.pcgtol,param.pcgmaxits); 
    end
            
end




end


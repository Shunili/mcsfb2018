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
%   reg_eps=1/G.N;
%     reg_eps=.001;
reg_eps=1;
else
    reg_eps=param.reg_eps;
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
    order=20;
else
    order=param.order;
end

if ~isfield(param,'grid_order')
    grid_order=1000;
else
    grid_order=param.grid_order;
end

%lower_1 = 0.5;
%upper_1 = 4.5;
lower_1=lower;
upper_1=upper;
h = @(x) (x>=lower_1 & x<upper_1);

if ~isfield(param,'reg_filter')
    param.reg_filter = 1;
else
    param.reg_filter = param.reg_filter;
end

% if param.reg_filter == 1
%     reg_filter =@(x) 1./(h(x)+reg_eps)-1/(1+reg_eps);
% else
%     reg_filter =@(x) 1-h(x);
% end
% reg_filter =@(x) 1./(h(x)+reg_eps)-1/(1+reg_eps);

wd = zeros(G.N,1);
wd(selected)=1./(G.N*weights);
%B=diag(wd);
B = spdiags(wd,0,G.N,G.N);
right_side = zeros(G.N,1);
right_side(selected)=values./(weights*G.N);


if (isfield(G,'U') && isfield(G,'e'))
    eig_inds=(h(G.e)~=0);
    LHS=G.U(selected,eig_inds);
    rec_coef=LHS\values;
    z=G.U(:,eig_inds)*rec_coef;
%    A=B+G.U*diag(reg_filter(G.e))*G.U'; 
%    z=A\right_side;
else
    
    switch param.reg_filter
        case 1 % rational 
            range=[0,G.lmax];
            grid_order=1000;
            
%             lower_wide=lower-(upper-lower)/4;
%             upper_wide=upper+(upper-lower)/4;
            
            [~, JCH]=gsp_jackson_cheby_coeff(lower, upper, range, order);
            reg_short=@(x) 1./(gsp_cheby_eval(x,JCH,range)+reg_eps);
            short_coeff=gsp_cheby_coeff(G,reg_short,200,grid_order);  % may need to use damping here
            %gte=@(x) gsp_cheby_op(G, JCH, x)+reg_eps*x;
            %LHS=@(z) B*z-(1/(1+reg_eps))*z+pcg(gte,z,1e-10,100);
            LHS=@(z) B*z-(1/(1+reg_eps))*z+gsp_cheby_op(G,short_coeff,z);
            %z=pcg(LHS,right_side,1e-10,100);
 
            %preconditioner=@(z) z./(wd+gamma);
            initial_guess=zeros(G.N,1);
            initial_guess(selected)=values;
            
            z=pcg(LHS,right_side,1e-10,100,[],[],initial_guess);     
        case 2 % 1-h
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
            %afun = @(x) B*x+gamma*gsp_filter(G,reg_filter,x,param);
            if precondition
                preconditioner=@(z) z./(wd+gamma);
                initial_guess=zeros(G.N,1);
                initial_guess(selected)=values;
                z=pcg(afun,right_side,param.pcgtol,param.pcgmaxits,preconditioner,[],initial_guess);
            else
                z=pcg(afun,right_side,param.pcgtol,param.pcgmaxits); 
            end
        otherwise % spline     
            reg_eps=1;
            k = 1/reg_eps;
            lower_wide=lower-(upper-lower)/4;
            upper_wide=upper+(upper-lower)/4;
            delta=.15*(upper-lower);
            num_per_side=4;
            tt=[0,(lower_wide-delta)/2,linspace(lower_wide-delta,lower_wide,num_per_side),linspace(upper_wide,upper_wide+delta,num_per_side),(G.lmax+upper_wide+delta)/2,G.lmax];
            if lower==0
                zero_pen=0;
            else
                zero_pen=((lower_wide-delta)/G.lmax+1)*k;
            end

            if upper>=G.lmax
                high_pen=0;
            else
                high_pen=(2-(upper_wide+delta)/G.lmax)*k;
            end

            ftt=[zero_pen,((lower_wide-delta)/(2*G.lmax)+1)*k,linspace(k,0,num_per_side),linspace(0,k,num_per_side),(1.5-(upper_wide+delta)/(2*G.lmax))*k,high_pen];
    %        xi=stieltjes_spline(tt,ftt); 
    %        LHS=@(z) B*z+stieltjes_op(G,z,xi,tt,100,1e-6); 
            pp = pchip(tt, ftt);
            pen=@(x)ppval(pp,x);
            penc=gsp_cheby_coeff(G,pen,order,grid_order);
            kk=1:order;
            damping_coeffs=((1-kk/(order+2))*sin(pi/(order+2)).*cos(kk*pi/(order+2))+(1/(order+2))*cos(pi/(order+2))*sin(kk*pi/(order+2)))/sin(pi/(order+2));
            damping_coeffs=[1,damping_coeffs]';
            penc=penc.*damping_coeffs;
            LHS=@(z) B*z+gsp_cheby_op(G,penc,z);

            initial_guess=zeros(G.N,1);
            initial_guess(selected)=values;
            preconditioner=@(z) z./(wd+1);
            z=pcg(LHS,right_side,1e-10,1000,preconditioner,[],initial_guess);
    end    
    
end




end


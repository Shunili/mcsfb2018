function [ z , JCH_reg ] = mcsfb_reconstruct_band2( G , selected, values, lower, upper, weights, param)
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

if ~isfield(param,'gamma')
    gamma=1; %1/G.N;
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
    order=30;
else
    order=param.order;
end

if ~isfield(param,'grid_order')
    grid_order=1000;
else
    grid_order=param.grid_order;
end

if ~isfield(param,'reg_filter')
    param.reg_filter = 3;
else
    param.reg_filter = param.reg_filter;
end

% if param.reg_filter == 1
%     reg_filter =@(x) 1./(h(x)+reg_eps)-1/(1+reg_eps);
% else
%     reg_filter =@(x) 1-h(x);
% end
% reg_filter =@(x) 1./(h(x)+reg_eps)-1/(1+reg_eps);

h = @(x) (x>=lower & x<upper);

wd = zeros(G.N,1);
wd(selected)=gamma./(weights);
B = spdiags(wd,0,G.N,G.N);
right_side = zeros(G.N,1);
right_side(selected)=gamma*values./(weights);

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
            if ~isfield(param,'reg_eps')
            %   reg_eps=1/G.N;
                reg_eps=(sqrt(5)-1)/2; % peak of penalty at 1/(eps(1+eps)=1 in this case
            else
                reg_eps=param.reg_eps;
            end
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
            preconditioner=@(z) z./(wd+1);
            initial_guess=zeros(G.N,1);
            initial_guess(selected)=values;
            z=pcg(LHS,right_side,param.pcgtol,param.pcgmaxits,preconditioner,[],initial_guess);     

        case 2 % 1-h
            range=[0,G.lmax];
            [~, JCH]=gsp_jackson_cheby_coeff(lower, upper, range, order);
            JCH_reg=zeros(size(JCH));
            JCH_reg(1)=2;
            JCH_reg=JCH_reg-JCH;
            reg_filter_eval= @(x)gsp_cheby_op(G,JCH_reg,x);
            afun = @(x) B*x+reg_filter_eval(x);
            preconditioner=@(z) z./(wd+1); % If preconditioner is a matrix, it should have similar diagonal elements as the LHS; if a function, it should apply M^{-1}; that is 1 over the diagonal of LHS
            initial_guess=zeros(G.N,1);
            initial_guess(selected)=values;
            z=pcg(afun,right_side,param.pcgtol,param.pcgmaxits,preconditioner,[],initial_guess);

        case 3 % spline     
            lower_wide=lower-(upper-lower)/4;
            upper_wide=upper+(upper-lower)/4;
            delta=.15*(upper-lower);
            num_per_side=4;
            tt=[0,(lower_wide-delta)/2,linspace(lower_wide-delta,lower_wide,num_per_side),linspace(upper_wide,upper_wide+delta,num_per_side),(G.lmax+upper_wide+delta)/2,G.lmax];
            if lower==0
                zero_pen=0;
            else
                zero_pen=((lower_wide-delta)/G.lmax+1);
            end

            if upper>=G.lmax
                high_pen=0;
            else
                high_pen=(2-(upper_wide+delta)/G.lmax);
            end

            ftt=[zero_pen,((lower_wide-delta)/(2*G.lmax)+1),linspace(1,0,num_per_side),linspace(0,1,num_per_side),(1.5-(upper_wide+delta)/(2*G.lmax)),high_pen];
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
            z=pcg(LHS,right_side,param.pcgtol,param.pcgmaxits,preconditioner,[],initial_guess);
            JCH_reg=penc;
        case 4 % approximate subspace by filtering random vectors 
            % does not currently work well. add a qr decomposition?
            % something else?
            range=[0,G.lmax];
            [~, JCH]=gsp_jackson_cheby_coeff(lower, upper, range, order);
            Filtered=gsp_cheby_opX(G,JCH);
            LHS=Filtered(selected,:);
            rec_coef=LHS\values;
            z=Filtered*rec_coef;
        otherwise
            error('Unknown reconstruction method');
    end    
    
end

end


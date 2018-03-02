function [ z ] = approx_reconstruct2( G , y, weights , selected, gamma , lower,upper , precondition)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 

h = @(x) (x>=lower & x<upper);
reg_filter =@(x) 1-h(x);
wd = zeros(G.N,1);
wd(selected)=1./weights(selected);
B=diag(wd);
right_side = zeros(G.N,1);
right_side(selected)=y(selected)./weights(selected);
 

if isfield(G,'U')
    A=B+gamma*G.U*diag(reg_filter(G.e))*G.U'; 
    z=A\right_side;
%     afun=@(x)A*x; % this is the only difference
%     eps_pc=.01;
%     param.grid_order=1000;
%     param.order=100;
%     pc_filter= @(x) 1./(gamma*(reg_filter(x)+eps_pc));
%     preconditioner=@(z) gsp_filter(G,pc_filter,z,param); 
%     z=pcg(afun,right_side,1e-14,5000,preconditioner);
else % use conjugate gradient to solve linear system of equations, need to be careful because Chebyshev polynomial approx can result in negative values, making A indefinite. Jackson coeffs should keep A positive definite
    %param.grid_order=1000;
    order=30;
    range=[0,G.lmax];
    [~, JCH]=gsp_jackson_cheby_coeff(lower, upper, range, order);
    JCH_reg=zeros(size(JCH));
    JCH_reg(1)=2;
    JCH_reg=JCH_reg-JCH;
    reg_filter_eval= @(x)gsp_cheby_op(G,JCH_reg,x);
    afun = @(x) B*x+gamma*reg_filter_eval(x);
%    afun = @(x) B*x+gamma*gsp_filter(G,reg_filter,x,param);
    if precondition
%         pc_filter_coeff=JCH;
%         pc_filter_coeff(1)=pc_filter_coeff(1)+2;
%         pc_filter_coeff=pc_filter_coeff/gamma;
%         pc_filter_eval=@(x)gsp_cheby_op(G,pc_filter_coeff,x);
%         z=pcg(afun,right_side,1e-14,5000,pc_filter_eval); 
        

        %eps_pc=.01;
        %pc_filter= @(x) 1./(gamma*(reg_filter_eval(x)+eps_pc));
        %preconditioner=@(z) gsp_filter(G,pc_filter,z,param); 
        %z=pcg(afun,right_side,1e-14,5000,preconditioner);    
        preconditioner=@(z) z./(wd+gamma);
        initial_guess=zeros(G.N,1);
        initial_guess(selected)=y(selected);
        z=pcg(afun,right_side,1e-14,5000,preconditioner,[],initial_guess);
    else
        z=pcg(afun,right_side,1e-14,5000); 
    end
        

end
 

 

end
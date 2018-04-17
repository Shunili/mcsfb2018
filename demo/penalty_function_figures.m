h = @(x) (x>=4 & x<6);
[~,filter_coeff]=gsp_jackson_cheby_coeff(4, 6,[0 10], 50);
approx_h=@(x) gsp_cheby_eval(x,filter_coeff,[0,10]);

reg_filters = cell(3,1);
reg_filters{1} = @(x) 1-h(x);
reg_filters{2} = @(x) (1-approx_h(x));

reg_eps = 0.01;
reg_filters{3} = @(x) 1./(approx_h(x)+reg_eps)-1/(1+reg_eps);

% figure; 
% hold on;
% xx = 0:0.001:10;
% plot(xx, h(xx),'LineWidth',3);
% plot(xx, approx_h(xx),'LineWidth',3);
% set(gca,'FontSize',24);

k = 1/reg_eps;
figure; 
hold on;
ylim([0,k]);
xx = 0:0.001:10;
plot(xx, k*reg_filters{1}(xx),'r:','LineWidth',3);
plot(xx, k*reg_filters{2}(xx),'LineWidth',3);
plot(xx, reg_filters{3}(xx),'g','LineWidth',3);
set(gca,'FontSize',24);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);
leg = legend('$\kappa(1-h_m(\lambda))$','$\kappa(1-\tilde{h}_m(\lambda))$','$\frac{1}{\tilde{h}_m(\lambda)+\epsilon}-\frac{1}{1+\epsilon}$');
set(leg,'Interpreter','latex');
set(leg,'FontSize',15);
set(leg,'Location','southwest');







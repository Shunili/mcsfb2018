clear all;
close all;

lmax=10;
lower=4;
upper=6;
order=50;
grid_order=1000;

h = @(x) (x>=lower & x<upper);
[~,filter_coeff]=gsp_jackson_cheby_coeff(lower, upper,[0 lmax], order);
approx_h=@(x) gsp_cheby_eval(x,filter_coeff,[0,lmax]);

reg_filters = cell(3,1);
reg_filters{1} = @(x) 1-h(x);
reg_filters{2} = @(x) (1-approx_h(x));

reg_eps = (sqrt(5)-1)/2;
reg_filters{3} = @(x) 1./(approx_h(x)+reg_eps)-1/(1+reg_eps);

% figure; 
% hold on;
% xx = 0:0.001:10;
% plot(xx, h(xx),'LineWidth',3);
% plot(xx, approx_h(xx),'LineWidth',3);
% set(gca,'FontSize',24);


% monotonic cubic spline
lower_wide=lower-.5;
upper_wide=upper+.5;
delta=.3;
num_per_side=4;
%tt=[(upper+lower)/2,linspace(upper,upper+delta,num_per_side),lmax];
%ftt=[0,linspace(0,k,num_per_side),(2-(upper+delta)/lmax)*k];


tt=[0,(lower_wide-delta)/2,linspace(lower_wide-delta,lower_wide,num_per_side),linspace(upper_wide,upper_wide+delta,num_per_side),(lmax+upper_wide+delta)/2,lmax];
ftt=[((lower_wide-delta)/lmax+1),((lower_wide-delta)/(2*lmax)+1),linspace(1,0,num_per_side),linspace(0,1,num_per_side),(1.5-(upper_wide+delta)/(2*lmax)),(2-(upper_wide+delta)/lmax)];
%[~,pp] = stieltjes_spline(tt,ftt);
pp = pchip(tt, ftt);
pen=@(x)ppval(pp,x);
penc=gsp_cheby_coeff([0,lmax],pen,order,grid_order);
kk=1:order;
damping_coeffs=((1-kk/(order+2))*sin(pi/(order+2)).*cos(kk*pi/(order+2))+(1/(order+2))*cos(pi/(order+2))*sin(kk*pi/(order+2)))/sin(pi/(order+2));
damping_coeffs=[1,damping_coeffs]';
penc=penc.*damping_coeffs;
spline_approx=@(x) gsp_cheby_eval(x,penc,[0,lmax]);            
            
figure; 
hold on;
ylim([0,1.5]);
xx = 0:0.001:10;
plot(xx, reg_filters{1}(xx),'r:','LineWidth',3);
plot(xx, reg_filters{2}(xx),'r','LineWidth',3);
plot(xx, reg_filters{3}(xx),'b','LineWidth',3);
plot(xx, pen(xx),'k:','LineWidth',3);
plot(xx, spline_approx(xx),'k','LineWidth',3);
set(gca,'FontSize',24);
xlabel('$\lambda$','Interpreter','LaTex','FontSize',24);
leg = legend('$1-h_m(\lambda)$','$1-\tilde{h}_m(\lambda)$','$\frac{1}{\tilde{h}_m(\lambda)+\epsilon}-\frac{1}{1+\epsilon}$','Cubic spline',['Poly approx of',char(10),'cubic spline']);
set(leg,'Interpreter','latex');
set(leg,'FontSize',15);
set(leg,'Location','southwest');







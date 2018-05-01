param.replacement = 0;
L=50;
[weights, P_min_half] = compute_sampling_weights(G,L,h_tilde);

% nb_meas = 50:50:500;
nb_meas=300:50:800;
num_trials=10;
num_m=length(nb_meas);
mean_squared_error = cell(3,1);
mean_squared_error{1} = zeros(num_m,1);
mean_squared_error{2} = zeros(num_m,1);
mean_squared_error{3} = zeros(num_m,1);

for i = 1:3
    for k = 1:num_m
        total_mse=0;
        for j = 1:num_trials
            [M, selected] = build_sampling_matrix(G, weights, nb_meas(k),param);
            analysis_coeffs = f(selected);
           
            synth_param.reg_filter = i;
            synth_param.reg_eps = 1e-2;
            synth_param.order = 80;
            f_reconstruct = mcsfb_reconstruct_band2(G, selected, analysis_coeffs, low_limit, up_limit, weights(selected), synth_param);
            error=abs(f-f_reconstruct);
            error2 = f-f_reconstruct;
            total_mse=total_mse+sum(error.^2)/G.N;
        end
        mean_squared_error{i}(k)=total_mse/num_trials;
    end
end

figure;
hold on;
plot(nb_meas,mean_squared_error{1},'LineWidth',4);
plot(nb_meas,mean_squared_error{2},'LineWidth',4);
plot(nb_meas,mean_squared_error{3},'LineWidth',4);
xlabel('number of samples', 'FontSize',24);
ylabel('Mean Squared Error', 'FontSize',24);
% xlim([50,500]);
% ylim([0,0.8]);

xlim([300,800]);
% ylim([0,0.8]);
set(gca,'FontSize',24);
set(gca,'box','off');
leg = legend('$(1-\tilde{h}_m(\lambda))$','$\frac{1}{\tilde{h}_m(\lambda)+\epsilon}-\frac{1}{1+\epsilon}$','Spline Reg Filter');
set(leg,'Interpreter','latex');
set(leg,'FontSize',15);
set(leg,'Location','northeast');

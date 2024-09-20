function plot_generated(dataset_name)
parameter_file = ['../GLM_data/', dataset_name,'/parameters.mat'];
load(parameter_file, "J", "h", "firing_rates");
load("par_sig.mat", "par_sig");
J_gen = sum(J, 3);
figure("Visible","off");

subplot(2,2,1)
plot(h);
title('h')

subplot(2,2,2)
plot(firing_rates{1});
title('relative firing rate')

subplot(2,2,3)
clim_gen = max(abs(J_gen), [], "all");
imagesc(J_gen);
clim([-clim_gen, clim_gen]);
colormap("jet");
colorbar;
axis square;

subplot(2,2,4)
N = size(par_sig, 1);
scatter(reshape(J_gen, [1, N*N]), reshape(par_sig(:, 2:end), [1, N*N]));


fig_path = ['../figures/GLM/', dataset_name];
check_path(fig_path);
fig_file = [fig_path, '/GLMparameters_' dataset_name, '_generated.png'];
saveas(gcf, fig_file);
end
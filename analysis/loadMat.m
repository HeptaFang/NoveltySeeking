
session_name = 'Simulated_mini_Mus_partial';
model_path = ['../GLM_model/', session_name, '/GLM_', session_name, '_1_DeltaPure_0_L2=2_2500.mat'];
load(model_path, 'model_par', 'N');
par_Mus = model_par;

session_name = 'Simulated_mini_Ctr_partial';
model_path = ['../GLM_model/', session_name, '/GLM_', session_name, '_1_DeltaPure_0_L2=2_2500.mat'];
load(model_path, 'model_par');
par_Ctr = model_par;

session_name = 'Simulated_mini_Syn_partial';
model_path = ['../GLM_model/', session_name, '/GLM_', session_name, '_1_DeltaPure_0_L2=2_2500.mat'];
load(model_path, 'model_par');
par_Syn = model_par;

f = figure();

for i=1:3
    subplot(3, 3, i*3-2);
    imagesc(par_Mus(:, ((i-1)*N+2):(i*N+1)));
    title(['Mus ', num2str(i)]);
    axis equal;
    clim([-1.5, 1.5]);

    subplot(3, 3, i*3-1);
    imagesc(par_Ctr(:, ((i-1)*N+2):(i*N+1)));
    title(['Ctr ', num2str(i)]);
    axis equal;
    clim([-1.5, 1.5]);

    subplot(3, 3, i*3);
    imagesc(par_Syn(:, ((i-1)*N+2):(i*N+1)));
    title(['Syn ', num2str(i)]);
    axis equal;
    clim([-1.5, 1.5]);
end
seed = 118;
rng(seed);

N = 30;
J = zeros(N, N, 3);
h = zeros(N, 1);
connection_density = 0.02;

% generate random J
for i = 1:N 
    for j = 1:N
        for k = 1:3
            if rand() < connection_density
                J(i, j, k) = randi([1, 4])*3 - 7.5;
            end
            % adjecent neurons have stronger connection
            if abs(i-j) <= 5 && rand() < 0.3
                J(i, j, k) = J(i, j, k) + 0.5 + rand();
            end
        end
    end
end

% generate random h (-6 to -4)
for i = 1:N
    h(i) = 2*rand() - 5;
end

h_higher = h + 0.5;
h_lower = h - 0.5;

n_trial = 400;
B = 2000;

% load kernels
load('../GLM_data/kernel_DeltaPure.mat', 'conn_kernels', 'kernel_len');

% generate random raster
rasters = cell(1, n_trial);
rasters_higher = cell(1, n_trial);
rasters_lower = cell(1, n_trial);
trial_len = ones(1, n_trial)*B;
spikes = cell(1, n_trial);
spikes_higher = cell(1, n_trial);
spikes_lower = cell(1, n_trial);
firing_rates = zeros(N, n_trial);
firing_rates_higher = zeros(N, n_trial);
firing_rates_lower = zeros(N, n_trial);

for trial_idx = 1:n_trial
    fprintf("trial %d\n", trial_idx);
    raster = zeros(N, B);
    raster_higher = zeros(N, B);
    raster_lower = zeros(N, B);
    for t = 1:B
        % vectorized version
        % copy h to htot
        htot = h;
        htot_higher = h_higher;
        htot_lower = h_lower;
        if t>=kernel_len
            for k = 1:3
                kernel = conn_kernels{k};
                raster_conv = raster(:, t:-1:t-kernel_len+1);
                raster_higher_conv = raster_higher(:, t:-1:t-kernel_len+1);
                raster_lower_conv = raster_lower(:, t:-1:t-kernel_len+1);
                htot = htot + reshape(sum(J(:, :, k)*(kernel.*raster_conv), 2), N, 1);
                htot_higher = htot_higher + reshape(sum(J(:, :, k)*(kernel.*raster_higher_conv), 2), N, 1);
                htot_lower = htot_lower + reshape(sum(J(:, :, k)*(kernel.*raster_lower_conv), 2), N, 1);
            end
        end
        lambda = exp(htot);
        lambda_higher = exp(htot_higher);
        lambda_lower = exp(htot_lower);
        p = lambda./(1+lambda);
        p_higher = lambda_higher./(1+lambda_higher);
        p_lower = lambda_lower./(1+lambda_lower);
        raster(:, t) = rand(N, 1) < p;
        raster_higher(:, t) = rand(N, 1) < p_higher;
        raster_lower(:, t) = rand(N, 1) < p_lower;
    end
    rasters{trial_idx} = raster;
    rasters_higher{trial_idx} = raster_higher;
    rasters_lower{trial_idx} = raster_lower;
    firing_rates(:, trial_idx) = sum(raster, 2)/B;
    firing_rates_higher(:, trial_idx) = sum(raster_higher, 2)/B;
    firing_rates_lower(:, trial_idx) = sum(raster_lower, 2)/B;
end

session_name_full = "Simulated";
cell_id = cell(1, N);
cuetype = zeros(1, n_trial);
channel = 1:N;
% borders = [-0.5];

check_path(['../GLM_data/', 'Simulated']);
check_path(['../GLM_data/', 'Simulated_higher']);
check_path(['../GLM_data/', 'Simulated_lower']);
save_file = ['../GLM_data/', 'Simulated/raster_Simulated_1_0.mat'];
save(save_file, 'N', 'rasters', 'J', 'h', "n_trial", "rasters", "session_name_full", "trial_len", "channel", "firing_rates");

session_name_full = "Simulated_higher";
save_file = ['../GLM_data/', 'Simulated_higher/raster_Simulated_higher_1_0.mat'];
h = h_higher;
rasters = rasters_higher;
firing_rates = firing_rates_higher;
save(save_file, 'N', 'rasters', 'J', 'h', "n_trial", "rasters", "session_name_full", "trial_len", "channel", "firing_rates");

session_name_full = "Simulated_lower";
save_file = ['../GLM_data/', 'Simulated_lower/raster_Simulated_lower_1_0.mat'];
h = h_lower;
rasters = rasters_lower;
firing_rates = firing_rates_lower;
save(save_file, 'N', 'rasters', 'J', 'h', "n_trial", "rasters", "session_name_full", "trial_len", "channel", "firing_rates");

borders = [10.5, 20.5];
save(['../GLM_data/Simulated/borders_Simulated_1.mat'], "borders");
save(['../GLM_data/Simulated_higher/borders_Simulated_higher_1.mat'], "borders");
save(['../GLM_data/Simulated_lower/borders_Simulated_lower_1.mat'], "borders");
function generate_GLM_raster(N, B, kernel_name)
%% Generate test raster
% load kernels
kernel_file = ['../GLM_data/', 'kernel_', kernel_name, '.mat'];
load(kernel_file, "conn_kernels", "PS_kernels", ...
    "n_conn_kernel", "n_PS_kernel", "kernel_len");


sparsity = 0.02;
base_h = -4.5;
noise = 0;

kernel_x = 1:30;
kernel = exp(-kernel_x/10);
max_delay = 30;

raster = zeros(N, B);
J = zeros(N, N, max_delay); % to, from, delay


for i = 1:N
    for j = 1:N
        if rand<sparsity && i~=j
            % random connection from j to i, -1 to 1.
            delay = randi(max_delay);
            amp = rand * 2 - 1;
            J(i, j, :) = amp*kernel;
        end
    end
end

%% manual parameters
h = base_h + rand(N, 1) - 0.5;
% high_fr_neurons
h(11:16) = -3;

% low_fr_neurons
h(31:33) = -7;


% global inhibitory neurons
for i = 1:30
    for j = 9:12
        if rand<0.3 && i~=j
            amp = rand - 1;
            J(i, j, :) = amp*kernel;
        end
    end
end

% global excitatory neurons
for i = 1:30
    for j = 15:18
        if rand<0.3 && i~=j
            amp = rand;
            J(i, j, :) = amp*kernel;
        end
    end
end

%% generate raster
h_tot = repmat(h, [1, B+max_delay]);
for t = max_delay+1:B
    % loop over time
    f_prob = 1-1./(1+exp(h_tot(:, t)));
    raster(:, t) = rand(N, 1) < f_prob;
    % delay connection
    h_tot(:, t+1:t+max_delay) = h_tot(:, t+1:t+max_delay) + ...
        tensorprod(J, raster(:, t), 2, 1);
end

% binary raster
raster(raster>1)=1;

firing_rate = mean(raster, 2);
firing_rates = {firing_rate};
rasters = {raster};
n_trial = 1;
trial_len = zeros(1, n_trial);
for i=1:n_trial
    trial_len(i) = size(rasters{i}, 2);
end

% save generated raster
dataset_name = 'generated';
session = 0;
check_path(['../GLM_data/', dataset_name]);
save(['../GLM_data/', dataset_name,'/parameters.mat'], "N", "B", "sparsity", "max_delay",...
    "base_h", "noise", "h", "J", "firing_rates");
save(['../GLM_data/', dataset_name,'/raster_',dataset_name, '_', int2str(session),'_0.mat'],...
    "rasters", "firing_rates", "n_trial", "trial_len");

end


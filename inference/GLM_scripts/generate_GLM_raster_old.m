% Generate test raster

% gaussian function
gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;

N = 50;
B = 44843;
sparsity = 0.02;
max_delay = 30;
delay_x = 1:max_delay;
base_h = -4.5;
noise = 0;

kernel_x = 1:30;
kernel = exp(-kernel_x/10);

raster = zeros(N, B);
J = zeros(N, N, max_delay); % to, from, delay
h = base_h + rand(N, 1) - 0.5;
h(1:5) = -3;
h_tot = repmat(h, [1, B+max_delay]);


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

for t = max_delay+1:B
    % loop over time
    f_prob = 1-1./(1+exp(h_tot(:, t)));
    raster(:, t) = rand(N, 1) < f_prob;
    % delay connection
    h_tot(:, t+1:t+max_delay) = h_tot(:, t+1:t+max_delay) + ...
        tensorprod(J, raster(:, t), 2, 1);
end

firing_rates = mean(raster, 2);
rasters = {raster};
n_trial = 1;
trial_len = zeros(n_trial, 1);
for i=1:n_trial
    trial_len(i) = size(rasters{i}, 2);
end
% save generated raster
dataset_name = 'generated';
session = 0;
check_path(['../GLM_data/', dataset_name]);
save(['../GLM_data/', dataset_name,'/parameters.mat'], "N", "B", "sparsity", "max_delay",...
    "base_h", "noise", "h", "J", "firing_rates");
save(['../GLM_data/', dataset_name,'/raster_',dataset_name, '_', int2str(session),'.mat'],...
    "rasters", "firing_rates", "n_trial", "trial_len");

% save("GLM_test_data/parameters.mat", "N", "B", "sparsity", "max_delay", "base_h", "noise",...
%     "h", "J", "firing_rates");
% save("GLM_test_data/raster.mat", "raster");
% Generate simulated data for the GLM model
% Validate possible model structures

%% Generate basic connectivity matrix
seed = 144; 
rng(seed);
% session_name_full = 'Simulated_Mus_partial';

% UseGPU = gpuDeviceCount > 0;
UseGPU = false;

groups = {'Thal', 'ACC-E', 'ACC-I', 'VLPFC-E', 'VLPFC-I'};
n_groups = length(groups);
n_cells = [10, 15, 5, 15, 5];
n_cells_cumulative = cumsum(n_cells);
n_cells_cumulative = [0, n_cells_cumulative];

subgroup_size = [15];
subgroup_size_cumulative = cumsum(subgroup_size);
subgroup_size_cumulative = [0, subgroup_size_cumulative];

global_J_multi = 2;
kernel_J_multi = [1, 0.5, 0.25];

% strength parameters
within_E = 2; 
across_E = 1;
TC_E = 1; 
TC_I = 2;
within_EI = 1; 
within_IE = -2;
across_EI = 0.5; 
across_IE = 0;

% baseline h parameters
thalamus_h = -3.5; % thalamus h
E_h = -4.5; % excitatory h
I_h = -4; % inhibitory h


% connection_density = [  0,      0,      0,      0,      0   ;
%                         0.3,    0.2,    0.15,    0.05,   0   ;
%                         0.6,    0.2,    0,      0.05,   0   ;
%                         0.3,    0.05,   0,      0.2,    0.15 ;
%                         0.6,    0.05,   0,      0.2,    0   ];
connection_density = [  0,      0,      0,      0,      0   ;
                        0.2,    0.3,    0.3,    0.05,   0   ;
                        0.4,    0.2,    0,      0.05,   0   ;
                        0.2,    0.05,   0,      0.3,    0.3 ;
                        0.8,    0.05,   0,      0.2,    0   ];
% connection_density = [  0,      0,      0,      0,      0   ;
%                         0,      0.2,    0.2,    0.05,   0   ;
%                         0,      0.2,    0,      0.05,   0   ;
%                         0,      0.05,   0,      0.2,    0.2 ;
%                         0,      0.05,   0,      0.2,    0   ];
connection_strength = [ 0,      0,      0,      0,      0   ;
                        TC_E,   within_E,   within_IE,  across_E,   0        ;
                        TC_I,   within_EI,  0,          across_EI,  0        ;
                        TC_E,   across_E,   0,          within_E,   within_IE;
                        TC_I,   across_EI,  0,          within_EI,  0        ];

baseline_h = [thalamus_h, E_h, I_h, E_h, I_h];

N = sum(n_cells);
J = zeros(N, N, 3);
h = zeros(N, 1);

% generate random J
for i = 1:N 
    area_i = find(n_cells_cumulative < i, 1, 'last');

    if area_i == 2 || area_i == 4
        subgroup_i = find(subgroup_size_cumulative < i-n_cells_cumulative(area_i), 1, 'last');
    else
        subgroup_i = 0;
    end

    for j = 1:N
        area_j = find(n_cells_cumulative < j, 1, 'last');

        if area_j == 2 || area_j == 4
            subgroup_j = find(subgroup_size_cumulative < j-n_cells_cumulative(area_j), 1, 'last');
        else
            subgroup_j = 0;
        end

        for k = 1:3
            if subgroup_i == subgroup_j && area_i == area_j && subgroup_i ~= 0
                % within subgroup connection
                density = connection_density(area_i, area_j) * 1;
                strength = connection_strength(area_i, area_j) * 1;
            else
                % between subgroup connection
                density = connection_density(area_i, area_j);
                strength = connection_strength(area_i, area_j);
            end

            if rand() < density
                J(i, j, k) = strength*(0.5 + rand());
            end
            % % adjecent neurons have stronger connection
            % if abs(i-j) <= 5 && rand() < 0.3
            %     J(i, j, k) = J(i, j, k) + 0.5 + rand();
            % end
        end
    end
end

J = J*global_J_multi;
for k = 1:3
    J(:, :, k) = J(:, :, k) - diag(diag(J(:, :, k))); % remove self connection
    J(:, :, k) = J(:, :, k) * kernel_J_multi(k); % scale by kernel_J_multi
end

% generate random h
for i = 1:N
    area_i = find(n_cells_cumulative < i, 1, 'last');
    h(i) = baseline_h(area_i) + randn() - 0.5;
end

n_trial = 50;
B = 10000;

% load kernels
load('../GLM_data/kernel_DeltaPure.mat', 'conn_kernels', 'kernel_len');

conn_kernels_reversed = cell(1, 3);
for k = 1:3
    conn_kernels_reversed{k} = conn_kernels{k}(:, end:-1:1);
end


%% generate random raster for control, muscimol and sync session

% session_names = {'Simulated_Ctr_partial_mini', 'Simulated_Sync_partial_mini', 'Simulated_Mus_partial_mini'};
session_names = {'Simulated_Ctr_mini', 'Simulated_Sync_partial_mini', 'Simulated_Mus_partial_mini'};
J_ctr = J;
all_firing_rates = cell(1, 3);
for session_idx = 1:3
    N = sum(n_cells);
    session_name_full = session_names{session_idx};

    if strcmp(session_name_full, 'Simulated_Mus_partial_mini')
        % remove thalamus connections
        J(:, 1:10, :) = 0;
        J_Mus = J;
    elseif strcmp(session_name_full, 'Simulated_Sync_partial_mini')
        % uniform thalamocortical connections
        sync_J = mean(J(11:50, 1:10, :), 2);
        J(11:50, 1:10, :) = repmat(sync_J, 1, 10, 1);
        J_Sync = J;
    end

    rasters = cell(1, n_trial);
    trial_len = ones(1, n_trial)*B;
    spikes = cell(1, n_trial);
    firing_rates = zeros(N, n_trial);

    % GPU if available
    if UseGPU
        fprintf("Using GPU\n");
        J = gpuArray(J);
        h = gpuArray(h);
        conn_kernels_reversed = cellfun(@gpuArray, conn_kernels_reversed, 'UniformOutput', false);
    end

    parfor trial_idx = 1:n_trial
        fprintf("%s, trial %d\n", session_name_full, trial_idx);
        raster = zeros(N, B);

        if UseGPU
            raster = gpuArray(raster);
        end

        for t = 1:B
            htot = h;
            if t>=kernel_len
                for k = 1:3
                    kernel = conn_kernels_reversed{k};
                    % raster_conv = raster(:, t-kernel_len+1:t);
                    % convolved = sum(raster_conv .* kernel, 2);
                    convolved = sum(raster(:, t-kernel_len+1:t) .* kernel, 2);
                    h_k = J(:, :, k)*convolved;
                    htot = htot + h_k;
                end
            end
            lambda = exp(htot);
            p = lambda./(1+lambda);
            raster(:, t) = rand(N, 1) < p;
        end

        if UseGPU
            raster = gather(raster);
        end

        rasters{trial_idx} = raster;
        firing_rates(:, trial_idx) = sum(raster, 2)/B;
    end

    all_firing_rates{session_idx} = firing_rates;

    % remove thalamus and inhibitory neurons
    for i=1:n_trial
        rasters{i}=[rasters{i}(11:25, :); rasters{i}(31:45, :)];
    end
    firing_rates = [firing_rates(11:25, :); firing_rates(31:45, :)];
    N = 30;

    cell_id = cell(1, N);
    cuetype = zeros(1, n_trial);
    channel = 1:N;
    % borders = [-0.5];

    check_path(['../GLM_data/', session_name_full]);
    save_file = ['../GLM_data/', session_name_full, '/raster_', session_name_full, '_1_0.mat'];
    save(save_file, 'N', 'rasters', 'J', 'h', "n_trial", "session_name_full", "trial_len", "channel", "firing_rates");

    % borders = [20.5, 60.5];
    borders = [15.5, 15.5];
    save(['../GLM_data/', session_name_full, '/borders_', session_name_full, '_1.mat'], "borders");

end

% plot J
all_J = {J_ctr, J_Sync, J_Mus};
figure('Position', [100, 100, 800, 800]);
t = tiledlayout(3, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
for i = 1:3
    for k = 1:3
        nexttile;
        imagesc(squeeze(all_J{i}(:, :, k)));
        % title(['J', num2str(k)]);
        colorbar;
        axis square;
        set(gca, 'XTick', [], 'YTick', []);
        title(['Session ', num2str(i), 'Kernel ', num2str(k)]);
    end
end

figure('Position', [100, 100, 800, 800]);

hold on;
for i = 1:3
    plot(mean(all_firing_rates{i}, 2), 'o-');
end
hold off;
title('Firing rates');
xlabel('Neuron ID');
ylabel('Firing rate (Hz)');
legend(session_names, 'Location', 'Best');
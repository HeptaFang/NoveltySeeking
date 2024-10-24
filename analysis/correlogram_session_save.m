function correlogram_session_save(session_type, session_idx, normalization)

% initialize session names
states = {...
    'PreDecision_full', 'PostDecision_cortex', 'PreInfo_full', 'PostInfo_cortex', ...
    'PreInfoAnti_full', 'PostInfoAnti_cortex','PreInfoResp_full', 'PostInfoResp_cortex',...
    'PreRestOpen_full', 'PostRestOpen_cortex','PreRestClose_full', 'PostRestClose_cortex',...
    };

dataset_names = cell(length(states), 1);
for i = 1:length(states)
    dataset_names{i} = [session_type, states{i}];
end

all_data = cell(length(states), 1);

session_cell_area_pre = NaN;
session_cell_id_pre = NaN;
session_N_pre = NaN;

session_cell_area_post = NaN;
session_cell_id_post = NaN;
session_N_post = NaN;

% load sessions
for i = 1:length(dataset_names)
    dataset_name = dataset_names{i};
    filename = ['../GLM_data/', dataset_name, '/raster_', dataset_name, '_', ...
        int2str(session_idx), '_0.mat'];
    load(filename, 'spikes'); % (N, n_trial) cells.
    load(filename, 'cell_area', 'cell_id', 'n_trial', 'N', 'trial_len');

    if strcmp(states{i}(1:3), 'Pre')
        if isnan(session_N_pre)
            session_cell_area_pre = cell_area;
            session_cell_id_pre = cell_id;
            session_N_pre = N;
        else
            assert(isequal(session_cell_area_pre, cell_area), 'cell_area not equal');
            assert(isequal(session_cell_id_pre, cell_id), 'cell_id not equal');
            assert(isequal(session_N_pre, N), 'N not equal');
        end
    else
        if isnan(session_N_post)
            session_cell_area_post = cell_area;
            session_cell_id_post = cell_id;
            session_N_post = N;
        else
            assert(isequal(session_cell_area_post, cell_area), 'cell_area not equal');
            assert(isequal(session_cell_id_post, cell_id), 'cell_id not equal');
            assert(isequal(session_N_post, N), 'N not equal');
        end
    end

    all_data{i}.spikes = spikes;
    all_data{i}.n_trial = n_trial;
    all_data{i}.trial_len = trial_len;
end

% check if pre and post session has the same cell id except thalamus
thalamus_idx_pre = strcmp(session_cell_area_pre, 'Thalamus');
thalamus_idx_post = strcmp(session_cell_area_post, 'Thalamus');
assert(isequal(session_cell_area_pre(~thalamus_idx_pre), ...
    session_cell_area_post(~thalamus_idx_post)), 'cell_area not equal');
assert(isequal(session_cell_id_pre(~thalamus_idx_pre), ...
    session_cell_id_post(~thalamus_idx_post)), 'cell_id not equal');

cell_area = session_cell_area_pre;
cell_id = session_cell_id_pre;
N = session_N_pre;
N_thalamus = sum(thalamus_idx_pre);

% parameters
temp_range = 200; % ms
N_shuffle = 100;

smooth_kernel = ones(1, 5) / 5; % rectangular kernel
% smooth_kernel = gausswin(5); % gaussian kernel

total_figure_num = N*(N-1)/2;

% for each dataset
for dataset_name_idx = 1:length(dataset_names)
    correlogram_map = zeros(N, N, temp_range+1);
    correlogram_map_shuffled = zeros(N, N, N_shuffle, temp_range+1);

    plot_count = 0;
    total_time_spent = 0;
    failed_count = 0;
    skipped_count = 0;
    failed_list = {};

    % for each neuron pair
    for i = 1:(N-1)
        for j = (i+1):N
            plot_count = plot_count + 1;
            fprintf("Calculating %d/%d\n", plot_count, total_figure_num);

            tic;

            try

                if (strcmp(cell_area{i}, 'Thalamus') || strcmp(cell_area{j}, 'Thalamus')) &&...
                    strcmp(states{dataset_name_idx}(1:4), 'Post')
                    continue;
                end

                if strcmp(cell_area{i}, 'VLPFC') && strcmp(states{dataset_name_idx}(1:4), 'Post')
                    i_fixed = i - N_thalamus;
                else
                    i_fixed = i;
                end

                if strcmp(cell_area{j}, 'VLPFC') && strcmp(states{dataset_name_idx}(1:4), 'Post')
                    j_fixed = j - N_thalamus;
                else
                    j_fixed = j;
                end

                spikes = all_data{dataset_name_idx}.spikes;
                n_trial = all_data{dataset_name_idx}.n_trial;
                trial_len = all_data{dataset_name_idx}.trial_len;

                dataset_name = dataset_names{dataset_name_idx};
                
                counts = zeros(1, 2*temp_range+1);
                counts_shuffled = zeros(N_shuffle, 2*temp_range+1);

                bin_edges = (-(temp_range+0.5)):1:(temp_range+0.5);
                bin_centers = (-temp_range):1:temp_range;
                total_len = sum(trial_len);
                total_spike_i = 0;
                total_spike_j = 0;
                chance_level = zeros(size(bin_centers));

                if strcmp(normalization, 'all')
                end
                for k = 1:n_trial
                    spikes1 = spikes{i_fixed, k};
                    spikes2 = spikes{j_fixed, k};    
                    
                    if isempty(spikes1) || isempty(spikes2)
                        continue;    
                    end
                    
                    switch normalization
                    case 'trial'
                        chance_trial = length(spikes1) * length(spikes2) *(trial_len(k) - abs(bin_centers)) / (trial_len(k)^2);
                        chance_level = chance_level + chance_trial;
                    case 'all'
                        total_spike_i = total_spike_i + length(spikes1);
                        total_spike_j = total_spike_j + length(spikes2);
                    end

                    trial_diffs = reshape(spikes1 - spikes2.', 1, []);
                    counts = counts + histcounts(trial_diffs*1000, bin_edges);
                    
                    shuffled_spikes1 = rand(N_shuffle, length(spikes1)) * trial_len(k)/1000;
                    shuffled_spikes2 = rand(N_shuffle, length(spikes2)) * trial_len(k)/1000;
                    for l = 1:N_shuffle
                        shuffled_diffs_trial = reshape(shuffled_spikes1(l, :) - shuffled_spikes2(l, :).', 1, []);
                        counts_shuffled(l, :) = counts_shuffled(l, :) + histcounts(shuffled_diffs_trial*1000, bin_edges);
                    end
                end

                % normalize counts by chance level
                switch normalization
                case 'trial'
                    % do nothing
                case 'all'
                    chance_level = total_spike_i * total_spike_j *(total_len - abs(bin_centers)) / (total_len^2);
                end
                count_over_chance = counts ./ chance_level;
                smoothed_count = conv(count_over_chance, smooth_kernel, 'same');
                
                % count each shuffled correlogram
                count_over_chance_shuffled = counts_shuffled ./ chance_level;
                smoothed_count_shuffled = zeros(N_shuffle, 2*temp_range+1);
                for k = 1:N_shuffle
                    smoothed_count_shuffled(k, :) = conv(count_over_chance_shuffled(k, :), smooth_kernel, 'same');
                end

                % log scale
                % count_over_chance_log = log(count_over_chance);
                % smoothed_count_log = log(smoothed_count);
                % mean_count_shuffled_log = log(mean_count_shuffled);
                % std_count_shuffled_log = std_count_shuffled/mean_count_shuffled;

                % save correlogram
                correlogram_map(i, j, :) = smoothed_count(temp_range+1:-1:1);
                correlogram_map(j, i, :) = smoothed_count(temp_range+1:end);
                correlogram_map_shuffled(i, j, :, :) = smoothed_count_shuffled(:, temp_range+1:-1:1);
                correlogram_map_shuffled(j, i, :, :) = smoothed_count_shuffled(:, temp_range+1:end);

            catch ME
                failed_count = failed_count + 1;
                failed_list{failed_count} = [cell_area{i}, '_', cell_id{i}, ...
                    '_vs_', cell_area{j}, '_', cell_id{j} ' reason: ', ME.message];
                fprintf("Failed: %s\n", failed_list{failed_count});
                close all;
            end

        fprintf("success: %d, skipped:%d, failed: %d\n", plot_count-skipped_count-failed_count, skipped_count, failed_count);

        time_spent = toc;
        total_time_spent = total_time_spent + time_spent;
        eta = (total_figure_num - plot_count)*total_time_spent/plot_count;
        fprintf("Time spent: %f, total: %02d:%02d:%02d, eta: %02d:%02d:%02d\n\n", time_spent, ...
            floor((total_time_spent/60/60)), ...
            floor(mod(total_time_spent/60, 60)), ...
            floor(mod(total_time_spent, 60)), ...
            floor((eta/60/60)), ...
            floor(mod(eta/60, 60)), ...
            floor(mod(eta, 60)));
        end
    end

    % save correlogram
    foldername = ['../GLM_data/', session_type, '_', int2str(session_idx)];
    check_path(foldername);

    model_path = [foldername, '/correlogramData_', dataset_name, '.mat'];

    save(model_path, 'correlogram_map', 'correlogram_map_shuffled', 'N', 'N_shuffle', ...
        'temp_range', 'cell_area', 'cell_id', 'normalization');
end

end
function correlogram_session(session_type, session_idx, normalization)

%%%%%%%%%%%%%%% todo: fix - post Muscimol sessions does not have thalamus data

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

% parameters
temp_range = 200; % ms
N_shuffle = 100;

smooth_kernel = ones(1, 5) / 5; % rectangular kernel
% smooth_kernel = gausswin(5); % gaussian kernel

total_figure_num = N*(N-1)/2;
plot_count = 0;
total_time_spent = 0;

% plotting for each neuron pair
for i = 1:(N-1)
    for j = (i+1):N
        plot_count = plot_count + 1;
        fprintf("Plot %d/%d\n", plot_count, total_figure_num);
        tic;
        % set up subplot 3x4
        figure("Visible", "off", "Position", [0, 0, 3200, 2400]);
        sgtitle([cell_area{i}, ' ', cell_id{i}, ' vs ', cell_area{j}, ' ', cell_id{j}]);
        tiles = tiledlayout(3, 4);

        global_max_y = 0;

        % loop over all sessions
        for dataset_name_idx = 1:length(dataset_names)
            nexttile;

            if (strcmp(cell_area{i}, 'Thalamus') || strcmp(cell_area{j}, 'Thalamus')) &&...
                strcmp(states{dataset_name_idx}(1:4), 'Post')
                continue;
            end

            spikes = all_data{dataset_name_idx}.spikes;
            n_trial = all_data{dataset_name_idx}.n_trial;
            trial_len = all_data{dataset_name_idx}.trial_len;

            dataset_name = dataset_names{dataset_name_idx};

            diffs = [];
            shuffled_diffs = zeros(N_shuffle, 0);
            bin_edges = (-(temp_range+0.5)):1:(temp_range+0.5);
            bin_centers = (-temp_range):1:temp_range;
            total_len = sum(trial_len);
            total_spike_i = 0;
            total_spike_j = 0;
            chance_level = zeros(size(bin_centers));

            if strcmp(normalization, 'all')
            end
            for k = 1:n_trial
                spikes1 = spikes{i, k};
                spikes2 = spikes{j, k};    
                
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


                diffs = [diffs, reshape(spikes1 - spikes2.', 1, [])];
                
                shuffled_spikes1 = rand(N_shuffle, length(spikes1)) * trial_len(k)/1000;
                shuffled_spikes2 = rand(N_shuffle, length(spikes2)) * trial_len(k)/1000;
                shuffled_diffs_trial = zeros(N_shuffle, length(spikes1)*length(spikes2));
                for l = 1:N_shuffle
                    shuffled_diffs_trial(l, :) = reshape(shuffled_spikes1(l, :) - shuffled_spikes2(l, :).', 1, []);
                end
                shuffled_diffs = [shuffled_diffs, shuffled_diffs_trial];
            end

            % Calculate histogram and plot curve
            counts = histcounts(diffs*1000, bin_edges);

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
            count_over_chance_shuffled = zeros(N_shuffle, length(bin_centers));
            smoothed_count_shuffled = zeros(N_shuffle, length(bin_centers));
            for k = 1:N_shuffle
                counts = histcounts(shuffled_diffs(k, :)*1000, bin_edges);
                count_over_chance_shuffled(k, :) = counts ./ chance_level;
                smoothed_count_shuffled(k, :) = conv(count_over_chance_shuffled(k, :), smooth_kernel, 'same');
            end

            % calculate mean and std of shuffled correlogram
            mean_count_shuffled = mean(count_over_chance_shuffled, 1);
            std_count_shuffled = std(count_over_chance_shuffled, 1);
            max_count_shuffled = max(count_over_chance_shuffled, [], 1);
            min_count_shuffled = min(count_over_chance_shuffled, [], 1);

            % log scale
            count_over_chance = log(count_over_chance);
            smoothed_count = log(smoothed_count);
            mean_count_shuffled = log(mean_count_shuffled);
            std_count_shuffled = std_count_shuffled/mean_count_shuffled;

            max_y = 1.1*max(max(abs(count_over_chance(~isinf(count_over_chance)))), 1);
            global_max_y = max(global_max_y, max_y);

            % plot correlogram
            hold on;
            % error bar area and mean
            fill([bin_centers, fliplr(bin_centers)], ...
                [mean_count_shuffled + std_count_shuffled, fliplr(mean_count_shuffled - std_count_shuffled)], ...
                [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
            % fill([bin_centers, fliplr(bin_centers)], ...
            %     [max_count_shuffled, fliplr(min_count_shuffled)], ...
            %     [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
            plot(bin_centers, mean_count_shuffled, 'k-', 'LineWidth', 1.5);

            % data
            plot(bin_centers, count_over_chance, 'r.-', 'LineWidth', 1.5, 'Color', [1, 0, 0, 0.5]);
            plot(bin_centers, smoothed_count, 'b-', 'LineWidth', 1.5);
            
            plot([-temp_range, temp_range], [0, 0], 'k--', 'LineWidth', 1);
            hold off;
            
            % replace underscore with space
            data_title = strrep(states{dataset_name_idx}, '_', ' ');

            title(data_title, 'FontSize', 18);
            xlabel('T_1 - T_2 (ms)');
            ylabel('log correlation');
            xlim([-temp_range, temp_range]);
        end

        % make all ylim the same
        for k = 1:12
            nexttile(k);
            ylim([-global_max_y, global_max_y]);
            hold on;
            plot([0, 0], [-global_max_y, global_max_y], 'k--', 'LineWidth', 1);
            hold off;
            legend('Shuffled range', 'Shuffled mean', 'Data', 'Smoothed data');
        end

        % global title
        sgtitle([cell_area{i}, ' ', cell_id{i}, ' vs ', cell_area{j}, ' ', cell_id{j}]);

        % save figure
        fig_folder = ['../figures/correlograms/', session_type, '_', int2str(session_idx)];
        if ~exist(fig_folder, 'dir')
            mkdir(fig_folder);
        end
        fig_name = [fig_folder, '/correlogram_', cell_area{i}, '_', cell_id{i}, ...
            '_vs_', cell_area{j}, '_', cell_id{j}, '.png'];
        fprintf("(%d/%d)Save %s\n", plot_count, total_figure_num, fig_name);
        saveas(gcf, fig_name);

        close all;
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

end
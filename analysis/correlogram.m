function correlogram(dataset_name, session, normalization)
% load data
filename = ['../GLM_data/', dataset_name, '/raster_', dataset_name, '_', ...
    int2str(session), '_0.mat'];
load(filename, 'spikes'); % (N, n_trial) cells.
load(filename, 'cell_area', 'cell_id', 'n_trial', 'N', 'trial_len');

temp_range = 200; % ms
N_shuffle = 100;

smooth_kernel = ones(1, 5) / 5; % rectangular kernel
% smooth_kernel = gausswin(5); % gaussian kernel

% calc spike time difference between all pairs of cells
for i = 1:N
    for j = (i+1):N
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

        % plot correlogram, bin size = 1ms, window = [-100ms, 100ms]
        figure("Visible", "off", "Position", [0, 0, 800, 600]);
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
        

        max_y = 1.1*max(max(count_over_chance), 1);
        hold on;
        % error bar area and mean
        % fill([bin_centers, fliplr(bin_centers)], ...
        %     [mean_count_shuffled + std_count_shuffled, fliplr(mean_count_shuffled - std_count_shuffled)], ...
        %     [0, 0, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        fill([bin_centers, fliplr(bin_centers)], ...
            [max_count_shuffled, fliplr(min_count_shuffled)], ...
            [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        plot(bin_centers, mean_count_shuffled, 'k-', 'LineWidth', 1.5);

        % data
        plot(bin_centers, count_over_chance, 'r.-', 'LineWidth', 1.5, 'Color', [1, 0, 0, 0.5]);
        plot(bin_centers, smoothed_count, 'b-', 'LineWidth', 1.5);
        
        plot([-temp_range, temp_range], [1, 1], 'k--', 'LineWidth', 1);
        plot([0, 0], [0, max_y], 'k--', 'LineWidth', 1);
        hold off;

        title([cell_area{i}, ' ', cell_id{i}, ' vs ', cell_area{j}, ' ', cell_id{j}]);
        xlabel('T_1 - T_2 (ms)');
        ylabel('Count over chance level');
        xlim([-temp_range, temp_range]);
        ylim([0, max_y]);
        legend('Shuffled range', 'Shuffled mean', 'Data', 'Smoothed data');

        % save figure
        fig_folder = ['../figures/correlograms/', dataset_name, '_', int2str(session)];
        if ~exist(fig_folder, 'dir')
            mkdir(fig_folder);
        end
        fig_name = [fig_folder, '/correlogram_', cell_area{i}, '_', cell_id{i}, ...
            '_vs_', cell_area{j}, '_', cell_id{j}, '.png'];
        fprintf("Save %s\n", fig_name);
        saveas(gcf, fig_name);

        close all;
    end
end
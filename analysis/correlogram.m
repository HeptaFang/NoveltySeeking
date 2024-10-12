function correlogram(dataset_name, session, normalization)
% load data
filename = ['../GLM_data/', dataset_name, '/raster_', dataset_name, '_', ...
    int2str(session), '_0.mat'];
load(filename, 'spikes'); % (N, n_trial) cells.
load(filename, 'cell_area', 'cell_id', 'n_trial', 'N', 'trial_len');


% calc spike time difference between all pairs of cells
for i = 1:N
    for j = (i+1):N
        diffs = [];
        bin_edges = -100.5:1:100.5;
        bin_centers = -100:1:100;
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
        
        max_y = 1.1*max(max(count_over_chance), 1);
        hold on;
        plot(bin_centers, count_over_chance, 'r.-', 'LineWidth', 1.5);
        plot([-100, 100], [1, 1], 'k--', 'LineWidth', 1);
        plot([0, 0], [0, max_y], 'k--', 'LineWidth', 1);
        hold off;

        title([cell_area{i}, ' ', cell_id{i}, ' vs ', cell_area{j}, ' ', cell_id{j}]);
        xlabel('T_1 - T_2 (ms)');
        ylabel('Count over chance level');
        xlim([-100, 100]);
        ylim([0, max_y]);

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
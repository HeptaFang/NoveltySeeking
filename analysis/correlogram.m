function correlogram(dataset_name, session)
% load data
filename = ['../GLM_data/', dataset_name, '/raster_', dataset_name, '_', ...
    int2str(session), '_0.mat'];
load(filename, 'spikes'); % (N, n_trial) cells.
load(filename, 'cell_area', 'cell_id', 'n_trial', 'N');

% calc spike time difference between all pairs of cells
for i = 1:N
    for j = (i+1):N
        diffs = [];
        for k = 1:n_trial
            spikes1 = spikes{i, k};
            spikes2 = spikes{j, k};
            if isempty(spikes1) || isempty(spikes2)
                continue;
            end
            diffs = [diffs, reshape(spikes1 - spikes2.', 1, [])];
        end

        % plot correlogram, bin size = 1ms, window = [-100ms, 100ms]
        figure;
        histogram(diffs*1000, -100:1:100);
        title([cell_area{i}, ' ', int2str(cell_id(i)), ' vs ', cell_area{j}, ' ', int2str(cell_id(j))]);
        xlabel('T_1 - T_2 (ms)');
        ylabel('Count');

        % save figure
        fig_folder = ['../figures/correlograms/', dataset_name, '_', int2str(session)];
        if ~exist(fig_folder, 'dir')
            mkdir(fig_folder);
        end
        fig_name = [fig_folder, '/correlogram_', cell_area{i}, '_', int2str(cell_id(i)), '_vs_', cell_area{j}, '_', int2str(cell_id(j)), '.png'];
        saveas(gcf, fig_name);

        close all;
    end
end
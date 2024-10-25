function correlogram_kernel_weighted_plot(session_type, session_idx)

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

%% load correlogram data
correlogram_map_all = NaN;
correlogram_map_shuffled_all = NaN;
correlogram_map_significance = NaN;
kernel_corr_all = NaN;
kernel_corr_significance_all = NaN;
N = NaN;
N_shuffle = NaN;
temp_range = NaN;
n_conn_kernel = NaN;

% Load data
fprintf('Loading data...\n');
for i = 1:length(states)
    dataset_name = dataset_names{i};
    foldername = ['../GLM_data/', session_type, '_', int2str(session_idx)];
    data_path = [foldername, '/correlogramData_', dataset_name, '.mat'];

    load_data = load(data_path, 'correlogram_map', 'correlogram_map_shuffled', 'correlogram_map_significance', ...
    'kernel_corr', 'kernel_corr_significance', 'conn_kernels', 'n_conn_kernel', 'N', 'N_shuffle', ...
    'temp_range', 'cell_area', 'cell_id', 'normalization');
    % data format:
    % correlogram_map = zeros(N, N, temp_range+1);
    % correlogram_map_shuffled = zeros(N, N, N_shuffle, temp_range+1);

    if isnan(N)
        N = load_data.N;
        N_shuffle = load_data.N_shuffle;
        temp_range = load_data.temp_range;
        n_conn_kernel = load_data.n_conn_kernel;
        correlogram_map_all = zeros(N, N, temp_range+1, length(states));
        correlogram_map_significance = zeros(N, N, temp_range+1, length(states));
        correlogram_map_shuffled_all = struct();
        correlogram_map_shuffled_all.mean = zeros(N, N, temp_range+1, length(states));
        correlogram_map_shuffled_all.std = zeros(N, N, temp_range+1, length(states));
        correlogram_map_shuffled_all.max = zeros(N, N, temp_range+1, length(states));
        correlogram_map_shuffled_all.min = zeros(N, N, temp_range+1, length(states));
        kernel_corr_all = zeros(N, N, n_conn_kernel, length(states));
        kernel_corr_significance_all = zeros(N, N, n_conn_kernel, length(states));
    else
        assert(N == load_data.N, 'N not equal');
        assert(N_shuffle == load_data.N_shuffle, 'N_shuffle not equal');
        assert(temp_range == load_data.temp_range, 'temp_range not equal');
        assert(n_conn_kernel == load_data.n_conn_kernel, 'n_conn_kernel not equal');
    end

    correlogram_map_all(:, :, :, i) = load_data.correlogram_map;
    correlogram_map_shuffled_all.mean(:, :, :, i) = load_data.correlogram_map_shuffled.mean;
    correlogram_map_shuffled_all.std(:, :, :, i) = load_data.correlogram_map_shuffled.std;
    correlogram_map_shuffled_all.max(:, :, :, i) = load_data.correlogram_map_shuffled.max;
    correlogram_map_shuffled_all.min(:, :, :, i) = load_data.correlogram_map_shuffled.min;
    correlogram_map_significance(:, :, :, i) = load_data.correlogram_map_significance;
    kernel_corr_all(:, :, :, i) = load_data.kernel_corr;
    kernel_corr_significance_all(:, :, :, i) = load_data.kernel_corr_significance;
end

%% Plotting
% grouped by kernel
for k = 1:n_conn_kernel
    fprintf('Plotting kernel %d/%d...\n', k, n_conn_kernel);
    figure('Visible', 'off', 'Position', [0, 0, 1600, 1200]);
    tiledlayout(3, 4);
    for i = 1:length(states)
        nexttile;

        % plot correlogram
        image_data = log10(kernel_corr_all(:, :, k, i));
        imagesc(image_data);
        colorbar;
        axis square;
        xlabel('From');
        ylabel('To');

        % generate title
        data_title = states{i};
        % replace underscore with space
        data_title = strrep(data_title, '_', ' ');
        if strcmp(data_title(1:3), 'Pre')
            data_title = ['Pre injection, ', data_title(4:end)];
        else
            data_title = ['Post injection, ', data_title(5:end)];
        end
        % replace short form with full form
        data_title = strrep(data_title, 'InfoResp', 'info cue to reward');
        data_title = strrep(data_title, 'InfoAnti', 'choice to info cue');
        data_title = strrep(data_title, 'Info', 'choice to reward');
        data_title = strrep(data_title, 'Decision', 'before choice');
        data_title = strrep(data_title, 'RestOpen', 'Resting state, eye-open');
        data_title = strrep(data_title, 'RestClose', 'Resting state, eye-close');
        % remove full and cortex
        data_title = strrep(data_title, ' full', '');
        data_title = strrep(data_title, ' cortex', '');
        title(data_title);

    end
    sgtitle(['Kernel ', num2str(k)]);

    file_path = ['../figures/kernel_weighted_correlogram/', session_type, '_', int2str(session_idx)];
    check_path(file_path);
    saveas(gcf, [file_path, '/kernel_', num2str(k), '.png']);
end

% grouped by state
for i = 1:length(states)
    fprintf('Plotting state %d/%d...\n', i, length(states));
    figure('Visible', 'off', 'Position', [0, 0, 1600, 1200]);
    tiledlayout(3, 4);
    for k = 1:n_conn_kernel
        nexttile;

        % plot correlogram
        image_data = log10(kernel_corr_all(:, :, k, i));
        imagesc(image_data);
        colorbar;
        axis square;
        xlabel('From');
        ylabel('To');

        title(['Kernel ', num2str(k)]);

    end
    data_title = states{i};
    % replace underscore with space
    data_title = strrep(data_title, '_', ' ');
    if strcmp(data_title(1:3), 'Pre')
        data_title = ['Pre injection, ', data_title(4:end)];
    else
        data_title = ['Post injection, ', data_title(5:end)];
    end
    % replace short form with full form
    data_title = strrep(data_title, 'InfoResp', 'info cue to reward');
    data_title = strrep(data_title, 'InfoAnti', 'choice to info cue');
    data_title = strrep(data_title, 'Info', 'choice to reward');
    data_title = strrep(data_title, 'Decision', 'before choice');
    data_title = strrep(data_title, 'RestOpen', 'Resting state, eye-open');
    data_title = strrep(data_title, 'RestClose', 'Resting state, eye-close');
    % remove full and cortex
    data_title = strrep(data_title, ' full', '');
    data_title = strrep(data_title, ' cortex', '');

    sgtitle(data_title);

    file_path = ['../figures/kernel_weighted_correlogram/', session_type, '_', int2str(session_idx)];
    check_path(file_path);
    file_title = strrep(data_title, ' ', '_');
    file_title = strrep(file_title, ',', '');

    saveas(gcf, [file_path, '/state_', file_title, '.png']);
end

end
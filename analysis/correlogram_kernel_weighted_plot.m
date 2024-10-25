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

% load correlogram data
correlogram_map = zeros(N, N, temp_range+1, length(states));
correlogram_map_shuffled = zeros(N, N, N_shuffle, temp_range+1, length(states));

for i = 1:length(states)
    % load data
    dataset_name = dataset_names{i};
    foldername = ['../GLM_data/', session_type, '_', int2str(session_idx)];
    data_path = [foldername, '/correlogramData_', dataset_name, '.mat'];

    load(data_path, 'correlogram_map', 'correlogram_map_shuffled', 'N', 'N_shuffle', ...
        'temp_range', 'cell_area', 'cell_id', 'normalization');
        % data format:
        % correlogram_map = zeros(N, N, temp_range+1);
        % correlogram_map_shuffled = zeros(N, N, N_shuffle, temp_range+1);
end
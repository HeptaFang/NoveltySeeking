data = load("../combined_task.mat").aligned_data;
dataset_names = {'eyeClose', 'eyeOpen'};

% extract rasters
for session_idx = 1:5
    n_trial = 1;
    N = data.n_cell(session_idx);
    predj = zeros(N, 0);
    raster = zeros(N, 0);
    trial_t = zeros(N, 0);
    rasters = cell(1, n_trial);
    firing_rates = cell(1, n_trial);
    trial_len = zeros(1, n_trial);
    subdatas = {data.spk_close, data.spk_open};

    for sub_idx = 1:2
        dataset_name = dataset_names{sub_idx};
        subdata = subdatas{sub_idx};
        N = sum(~cellfun(@isempty, subdata(:, session_idx)));

        % fprintf('session %d, trial %d/%d\n', session_idx, trial_idx, n_trial);
        B = ceil(max(cellfun(@max, subdata(1:N, session_idx))));
        trial_raster = zeros(N, B);

        % calc raster and predj of this trial
        edges = 0:B;
        for i=1:N
            trial_raster(i, :) = histcounts(subdata{i, session_idx},edges);
        end

        trial_raster(trial_raster>1) = 1;

        rasters{1} = trial_raster;
        trial_len(1) = B;
        firing_rates{1} = mean(trial_raster, 2);

        check_path(['../GLM_data/', dataset_name]);
        save(['../GLM_data/', dataset_name,'/raster_', dataset_name, '_', int2str(session_idx),'_0.mat'],...
        "rasters", "firing_rates", "n_trial", "trial_len");
    end
    
end

% extract area borders
sub_idxes = {data.spk_close_idx, data.spk_open_idx};
for sub_idx = 1:2
    dataset_name = dataset_names{sub_idx};
    subdata = subdatas{sub_idx};
    neuron_idx = sub_idxes{sub_idx};
    neuron_idx = cellfun(@mean, neuron_idx);
    for session_idx=1:5
        borders = zeros(1,2);
        A_T_border = find(neuron_idx(:, session_idx)>45, 1);
        T_P_border = find(neuron_idx(:, session_idx)>98, 1);
        if isempty(A_T_border)
            A_T_border = max_n + 1;
        end
        if isempty(T_P_border)
            T_P_border = max_n + 1;
        end
        A_T_border = A_T_border-0.5;
        T_P_border = T_P_border-0.5;

        borders(1) = A_T_border;
        borders(2) = T_P_border;
        check_path(['../GLM_data/', dataset_name]);
        save(['../GLM_data/', dataset_name,'/borders_', dataset_name, '_', ...
            int2str(session_idx),'.mat'], "borders");
    end
end

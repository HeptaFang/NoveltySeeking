data = load("../combined_task.mat").aligned_data;
dataset_name = 'taskCue';

cue_max_t = 1009;
iti_max_t = 764;

% extract rasters
for session_idx = 1:5
    n_trial_all = length(data.cuetype{session_idx});
    N = data.n_cell(session_idx);
    predj = zeros(N, 0);
    raster = zeros(N, 0);
    trial_t = zeros(N, 0);
    rasters_all = cell(1, n_trial_all);
    firing_rates_all = cell(1, n_trial_all);
    trial_len_all = zeros(1, n_trial_all);

    for trial_idx = 1:n_trial_all
        fprintf('session %d, trial %d/%d\n', session_idx, trial_idx, n_trial_all);
        B = cue_max_t;
        trial_raster = zeros(N, B);

        % calc raster and predj of this trial
        edges = 0:0.001:B*0.001;
        for i=1:N
            trial_raster(i, :) = histcounts(data.spk_cue{i, session_idx}{trial_idx},edges);
        end

        trial_raster(trial_raster>1) = 1;

        rasters_all{trial_idx} = trial_raster;
        trial_len_all(trial_idx) = B;
        firing_rates_all{trial_idx} = mean(trial_raster, 2);
    end

    for cuetype = 1:5
        mask = data.cuetype{session_idx} == 6700+10*cuetype;
        rasters = rasters_all(mask);
        firing_rates = firing_rates_all(mask);
        n_trial = sum(mask);
        trial_len = trial_len_all(mask);
        session_idx_cue = session_idx*100+cuetype;

        check_path(['../GLM_data/', dataset_name]);
        save(['../GLM_data/', dataset_name,'/raster_', dataset_name, '_', int2str(session_idx_cue),'_0.mat'],...
        "rasters", "firing_rates", "n_trial", "trial_len");
    end
end

% extract area borders
area = data.area;
for session_idx=1:5
    borders = zeros(1,2);
    A_T_border = find(area(:, session_idx)==2, 1);
    T_P_border = find(area(:, session_idx)==3, 1);
    if isempty(A_T_border)
        A_T_border = N + 1;
    end
    if isempty(T_P_border)
        T_P_border = N + 1;
    end
    A_T_border = A_T_border-0.5;
    T_P_border = T_P_border-0.5;
    borders(1) = A_T_border;
    borders(2) = T_P_border;
    check_path(['../GLM_data/', dataset_name]);
    save(['../GLM_data/', dataset_name,'/borders_', dataset_name, '_', ...
        int2str(session_idx),'.mat'], "borders");
end


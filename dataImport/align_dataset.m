% Align task, eyeclose and eyeopen datasets, keep data length (B) the same.

%% Load data
control = 'Muscimol';
kernel = 'DeltaPure';

tasks = {'PreTask_full', 'PreTask_cortex', 'PostTask_cortex',...
    'PreRestClose_full', 'PreRestClose_cortex', 'PostRestClose_cortex',...
    'PreRestOpen_full', 'PreRestOpen_cortex', 'PostRestOpen_cortex'};

% check data size
for session = 1:10
    for task_idx = 1:length(tasks)
        task = tasks{task_idx};

        file_path = ['../GLM_data/', control, task, '/GLMdata_', control, task, '_', int2str(session), '_', kernel, '_0.mat'];
        load(file_path, 'N', 'B');
        file_path = ['../GLM_data/', control, task, '/raster_', control, task, '_', int2str(session), '_0.mat'];
        load(file_path, 'n_trial');
        
        fprintf('Session %d, task %s, N = %d, B = %d, n_trial = %d\n', session, task, N, B, n_trial);
    end
end

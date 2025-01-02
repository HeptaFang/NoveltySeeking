% Align task, eyeclose and eyeopen datasets, keep data length (B) the same.

%% Load data
control = 'Muscimol';
kernel = 'DeltaPure';

tasks = {'PreTask_full', 'PreTask_cortex', 'PostTask_cortex',...
    'PreRestClose_full', 'PreRestClose_cortex', 'PostRestClose_cortex',...
    'PreRestOpen_full', 'PreRestOpen_cortex', 'PostRestOpen_cortex'};

% check data size
for session = 1:10
    min_duration = 10000000;
    for task_idx = 1:length(tasks)
        task = tasks{task_idx};

        file_path = ['../GLM_data/', control, task, '/GLMdata_', control, task, '_', int2str(session), '_', kernel, '_0.mat'];
        load(file_path, 'N', 'B');
        file_path = ['../GLM_data/', control, task, '/raster_', control, task, '_', int2str(session), '_0.mat'];
        load(file_path, 'n_trial');

        min_duration = min(min_duration, B);
        
        fprintf('Session %d, task %s, N = %d, B = %d, n_trial = %d\n', session, task, N, B, n_trial);
    end
    seconds = floor(min_duration/1000);
    fprintf('Session %d, min_duration = %d, %d:%d\n', session, min_duration, floor(seconds/60), mod(seconds, 60));

    % align data to min_duration
    modes = {'First', 'Last', 'Random'};
    for task_idx = 1:length(tasks)
        task = tasks{task_idx};
        fprintf('Aligning %s %s\n', control, task);
        fprintf('Loading...');

        % raster file and border file
        file_path = ['../GLM_data/', control, task, '/raster_', control, task, '_', int2str(session), '_0.mat'];
        load(file_path, 'n_trial', 'cell_id', 'cell_area', 'channel', 'session_name_full');
        file_path = ['../GLM_data/', control, task, '/borders_', control, task, '_', int2str(session), '.mat'];
        load(file_path, 'borders');
        fprintf('Done\n');

        for mode_idx = 1:length(modes)
            mode = modes{mode_idx};
            fprintf('Mode: %s\n', mode);

            full_name = [control, task, '_Align', mode];
            fprintf('Saving...');
            check_path(['../GLM_data/', full_name]);

            % raster file and border file
            file_path = ['../GLM_data/', full_name, '/raster_', full_name, '_', int2str(session), '_0.mat'];
            save(file_path, 'n_trial', 'cell_id', 'cell_area', 'channel', 'session_name_full');
            file_path = ['../GLM_data/', full_name, '/borders_', full_name, '_', int2str(session), '.mat'];
            save(file_path, 'borders');

            fprintf('Done.\n');
        end
    end
end

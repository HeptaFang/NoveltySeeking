clear

task_names = {'Muscimol', 'Saline', 'SimRec'};

force_retrain = false;
total_training = 0;
skipped = 0;
failed= 0;
success = 0;

failed_list = {};

for task_idx=1:3
    fprintf("Task %d\n", task_idx);
    % for cuetype=1:5
    % compare if is Muscimol sessions
    if strcmp(task_names{task_idx}, 'Muscimol')
        session_idxs = 1:10;
    else
        session_idxs = 1:5;
    end

    for session_idx = session_idxs
        
        % try
            tick_session = tic;
            dataset_name = task_names{task_idx};
            fprintf("Main: %s, session%d\n", dataset_name, session_idx);
            total_training = total_training + 1;
            skip_flag = true;
            
            tic;
            correlogram_session_save(dataset_name, session_idx, 'trial');
            toc;

            success = success + 1;
        
        % catch ME
        %     fprintf("Failed: %s\n", ME.message);
        %     failed = failed + 1;
        %     failed_list{end+1} = {dataset_name, int2str(session_idx), ME.message};
        %     throw(ME);
        % end
    end
    fprintf("Task %d done\n", task_idx);
end

fprintf("Total: %d, Success: %d, Skipped: %d, Failed: %d\n", total_training, success, skipped, failed);
% save failed_list
save('../GLM_data/failed_list_corr.mat', 'failed_list');
if failed>0
    for i=1:length(failed_list)
        fprintf("Failed: %s, %s, %s\n", failed_list{i}{1}, failed_list{i}{2}, failed_list{i}{3});
    end
end
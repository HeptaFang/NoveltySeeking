%% GLM raster fitting
% Input: raw raster file "../GLM_data/[dataset_name]/raster_[dataset_name]_[session]_0.mat"

% Required variables in the input file:
% "n_trial": integer, trial number.
% "trial_len": int(1, n_trial), time bin number of each trial.
% "rasters": cell(1, n_trial), each element is a trial.
% Each raster is a (N, trial_len(i)) binary matrix, N is number of neurons.

clear
% task_names = {'MuscimolPre_cortex', 'MuscimolPost_cortex', 'MuscimolPre_full'};
% task_names = {...
%     'MuscimolPreDecision_cortex', 'MuscimolPostDecision_cortex', 'MuscimolPreInfo_cortex', 'MuscimolPostInfo_cortex', ...
%     'MuscimolPreInfoAnti_cortex', 'MuscimolPostInfoAnti_cortex','MuscimolPreInfoResp_cortex', 'MuscimolPostInfoResp_cortex',...
%     'MuscimolPreRestOpen_cortex', 'MuscimolPostRestOpen_cortex','MuscimolPreRestClose_cortex', 'MuscimolPostRestClose_cortex',...
%     'SalinePreDecision_cortex', 'SalinePostDecision_cortex', 'SalinePreInfo_cortex', 'SalinePostInfo_cortex', ...
%     'SalinePreInfoAnti_cortex', 'SalinePostInfoAnti_cortex','SalinePreInfoResp_cortex', 'SalinePostInfoResp_cortex',...
%     'SalinePreRestOpen_cortex', 'SalinePostRestOpen_cortex','SalinePreRestClose_cortex', 'SalinePostRestClose_cortex',...
%     };
task_names = {...
    'MuscimolPreDecision_full', 'MuscimolPostDecision_cortex', 'MuscimolPreInfo_full', 'MuscimolPostInfo_cortex', ...
    'MuscimolPreInfoAnti_full', 'MuscimolPostInfoAnti_cortex','MuscimolPreInfoResp_full', 'MuscimolPostInfoResp_cortex',...
    'MuscimolPreRestOpen_full', 'MuscimolPostRestOpen_cortex','MuscimolPreRestClose_full', 'MuscimolPostRestClose_cortex',...
    'SalinePreDecision_full', 'SalinePostDecision_cortex', 'SalinePreInfo_full', 'SalinePostInfo_cortex', ...
    'SalinePreInfoAnti_full', 'SalinePostInfoAnti_cortex','SalinePreInfoResp_full', 'SalinePostInfoResp_cortex',...
    'SalinePreRestOpen_full', 'SalinePostRestOpen_cortex','SalinePreRestClose_full', 'SalinePostRestClose_cortex',...
    'SimRecPreDecision_full', 'SimRecPreInfo_full',...
    'SimRecPreInfoAnti_full', 'SimRecPreInfoResp_full',...
    'SimRecPreRestOpen_full', 'SimRecPreRestClose_full',...
    };
force_retrain = false;
total_training = 0;
skipped = 0;
failed= 0;
success = 0;

failed_list = {};

% task_names = {'MuscimolPreDecision_full', 'SalinePreDecision_full', 'SimRecPreDecision_full'};
% task_names = {'MuscimolPostDecision_full', 'SalinePostDecision_full'};
% trial_names = {'100B', '50BI', '50BN', '100S', '0'};
for task_idx=19:20
    fprintf("Task %d\n", task_idx);
    % for cuetype=1:5
    % compare if is Muscimol sessions
    if strcmp(task_names{task_idx}(1:8), 'Muscimol')
        session_idxs = [6,7,8,9,10,1,4,5,2,3];
    else
        session_idxs = [1,4,5,2,3];
    end

    for session_idx = session_idxs
        for trial_idx = 1:1
            % try
                % if task_idx==1 && session_idx<10
                %     continue;
                % end
                % session_idx = session*100+cuetype;
                % fprintf("Main: session%d, cue%d\n", session, cuetype);

                % dataset_name = [task_names{task_idx}, '_', trial_names{trial_idx}];
                tick_session = tic;
                dataset_name = task_names{task_idx};
                fprintf("Main: %s, session%d\n", dataset_name, session_idx);
                total_training = total_training + 1;
                skip_flag = true;

                %% parameters
                % dataset_name = 'generated';
                % % session = 0;
                % kernel_name = 'expDecay10';
                
                % dataset_name = 'taskCue';
                % % session = 1;
                % kernel_name = 'expDecay10';
                
                tic;
                correlogram(dataset_name, session_idx, 'trial');
                toc;

                success = success + 1;
            
            % catch ME
            %     fprintf("Failed: %s\n", ME.message);
            %     failed = failed + 1;
            %     failed_list{end+1} = {dataset_name, int2str(session_idx), ME.message};
            %     throw(ME);
            % end
        end
    end
    % plot_rest(session_idx, kernel_name, max_epoch, reg, shuffle_size);
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
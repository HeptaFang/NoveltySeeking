%% GLM raster fitting
% Input: raw raster file "../GLM_data/[dataset_name]/raster_[dataset_name]_[session]_0.mat"

% Required variables in the input file:
% "n_trial": integer, trial number.
% "trial_len": int(1, n_trial), time bin number of each trial.
% "rasters": cell(1, n_trial), each element is a trial.
% Each raster is a (N, trial_len(i)) binary matrix, N is number of neurons.

clear
% task_names = {'MuscimolPre_cortex', 'MuscimolPost_cortex', 'MuscimolPre_full'};
task_names = {...
    % 'Simulated', 'Simulated_higher', 'Simulated_lower', ...
    % 'MuscimolPreOffer1_cortex', 'MuscimolPostOffer1_cortex', 'MuscimolPreOffer2_cortex', 'MuscimolPostOffer2_cortex', ...
    % 'MuscimolPreDecision_cortex', 'MuscimolPostDecision_cortex', 'MuscimolPreInfoAnti_cortex', 'MuscimolPostInfoAnti_cortex', ...
    % 'MuscimolPreInfoResp_cortex', 'MuscimolPostInfoResp_cortex', 'MuscimolPreReward_cortex', 'MuscimolPostReward_cortex', ...
    % 'MuscimolPreRandomA_cortex', 'MuscimolPostRandomA_cortex', 'MuscimolPreRandomB_cortex', 'MuscimolPostRandomB_cortex', ...
    % 'MuscimolPreOffer1_full',  'MuscimolPreOffer2_full',...
    % 'MuscimolPreDecision_full', 'MuscimolPreInfoAnti_full',...
    % 'MuscimolPreInfoResp_full', 'MuscimolPreReward_full',...
    % 'MuscimolPreRandomA_full',  'MuscimolPreRandomB_full',...
    % 'MuscimolPreRestOpen_cortex',...
    % 'MuscimolPostRestOpen_cortex', ...
    % 'MuscimolPreRestClose_cortex', ...
    % 'MuscimolPostRestClose_cortex', ...
    % 'MuscimolPreRestOpen_full',...
    % 'MuscimolPreRestClose_full', ...
    % 'MuscimolPreTask_full', 'MuscimolPostTask_cortex','MuscimolPreTask_cortex', ...

    % 'MuscimolPostTask_cortex_AlignRandom','MuscimolPreTask_cortex_AlignRandom', ...
    % 'MuscimolPostRestOpen_cortex_AlignRandom', 'MuscimolPreRestOpen_cortex_AlignRandom', ...
    % 'MuscimolPostRestClose_cortex_AlignRandom', 'MuscimolPreRestClose_cortex_AlignRandom', ...
    % 'MuscimolPostTask_cortex_AlignFirst','MuscimolPreTask_cortex_AlignFirst', ...
    % 'MuscimolPreTask_full_AlignRandom', 'MuscimolPreRestOpen_full_AlignRandom', 'MuscimolPreRestClose_full_AlignRandom', ...
    % 'MuscimolPreTask_full_AlignFirst', 'MuscimolPreRestOpen_full_AlignFirst', 'MuscimolPreRestClose_full_AlignFirst', ...
    % 
    % 'MuscimolPostTask_cortex_AlignRandom','MuscimolPreTask_cortex_AlignRandom', ...
    % 'MuscimolPostRestOpen_cortex_AlignRandom', 'MuscimolPreRestOpen_cortex_AlignRandom', ...
    % 'MuscimolPostRestClose_cortex_AlignRandom', 'MuscimolPreRestClose_cortex_AlignRandom', ...
    % 'MuscimolPostTask_cortex_AlignFirst','MuscimolPreTask_cortex_AlignFirst', ...
    % 'MuscimolPreTask_full_AlignRandom', 'MuscimolPreRestOpen_full_AlignRandom', 'MuscimolPreRestClose_full_AlignRandom', ...
    % 'MuscimolPreTask_full_AlignFirst', 'MuscimolPreRestOpen_full_AlignFirst', 'MuscimolPreRestClose_full_AlignFirst', ...
    % 'MuscimolPostRestOpen_cortex_AlignFirst', 'MuscimolPreRestOpen_cortex_AlignFirst', ...
    % 'MuscimolPostRestClose_cortex_AlignFirst', 'MuscimolPreRestClose_cortex_AlignFirst', ...

    % 'MuscimolPostRestOpen_cortex_AlignLast', 'MuscimolPreRestOpen_cortex_AlignLast', ...
    % 'MuscimolPostRestClose_cortex_AlignLast', 'MuscimolPreRestClose_cortex_AlignLast', ...
    % 'MuscimolPreRestOpen_full_AlignLast', 'MuscimolPreRestClose_full_AlignLast', ...

    % 'SalinePostRestOpen_cortex_AlignLast', 'SalinePreRestOpen_cortex_AlignLast', ...
    % 'SalinePostRestClose_cortex_AlignLast', 'SalinePreRestClose_cortex_AlignLast', ...
    % 'SalinePreRestOpen_full_AlignLast', 'SalinePreRestClose_full_AlignLast', ...
    'MuscimolPreAll_cortex_AlignLast', 'MuscimolPostAll_cortex_AlignLast','MuscimolPreAll_full_AlignLast',...
    'SalinePreAll_cortex_AlignLast', 'SalinePostAll_cortex_AlignLast','SalinePreAll_full_AlignLast',...
    
    % 'MuscimolPostTask_cortex_AlignLast','MuscimolPreTask_cortex_AlignLast', 'MuscimolPreTask_full_AlignLast', ...

    % 'SalinePostTask_cortex_AlignRandom','SalinePreTask_cortex_AlignRandom', ...
    % 'SalinePostRestOpen_cortex_AlignRandom', 'SalinePreRestOpen_cortex_AlignRandom', ...
    % 'SalinePostRestClose_cortex_AlignRandom', 'SalinePreRestClose_cortex_AlignRandom', ...
    % 'SalinePostTask_cortex_AlignFirst','SalinePreTask_cortex_AlignFirst', ...
    % 'SalinePreTask_full_AlignRandom', 'SalinePreRestOpen_full_AlignRandom', 'SalinePreRestClose_full_AlignRandom', ...
    % 'SalinePreTask_full_AlignFirst', 'SalinePreRestOpen_full_AlignFirst', 'SalinePreRestClose_full_AlignFirst', ...

    % 'MuscimolPostRestOpen_cortex_AlignFirst', 'MuscimolPreRestOpen_cortex_AlignFirst', ...
    % 'MuscimolPostRestClose_cortex_AlignFirst', 'MuscimolPreRestClose_cortex_AlignFirst', ...
    % 'MuscimolPostTask_cortex_AlignLast','MuscimolPreTask_cortex_AlignLast', ...
    % 'MuscimolPostRestOpen_cortex_AlignLast', 'MuscimolPreRestOpen_cortex_AlignLast', ...
    % 'MuscimolPostRestClose_cortex_AlignLast', 'MuscimolPreRestClose_cortex_AlignLast', ...
    % 'MuscimolPreTask_full_AlignLast', 'MuscimolPreRestOpen_full_AlignLast', 'MuscimolPreRestClose_full_AlignLast', ...
    % 'MuscimolPreRandomShort_full',  'MuscimolPreRandomLong_full',...
    % 'MuscimolPostRandomShort_cortex', 'MuscimolPostRandomLong_cortex',...
    % 'MuscimolPreRandomShort_cortex', 'MuscimolPreRandomLong_cortex', ...
    % 'SalinePreDecision_cortex', 'SalinePostDecision_cortex', 'SalinePreInfo_cortex', 'SalinePostInfo_cortex', ...
    % 'SalinePreInfoAnti_cortex', 'SalinePostInfoAnti_cortex','SalinePreInfoResp_cortex', 'SalinePostInfoResp_cortex',...
    % 'SalinePreRestOpen_cortex', 'SalinePostRestOpen_cortex','SalinePreRestClose_cortex', 'SalinePostRestClose_cortex',...
    };
    % task_names = {...
    %     'MuscimolPreDecision_full', 'MuscimolPostDecision_cortex', 'MuscimolPreInfo_full', 'MuscimolPostInfo_cortex', ...
    %     'MuscimolPreInfoAnti_full', 'MuscimolPostInfoAnti_cortex','MuscimolPreInfoResp_full', 'MuscimolPostInfoResp_cortex',...
    %     'MuscimolPreRestOpen_full', 'MuscimolPostRestOpen_cortex','MuscimolPreRestClose_full', 'MuscimolPostRestClose_cortex',...
    %     'MuscimolPreRestOpen_full', 'MuscimolPostRestOpen_cortex','MuscimolPreRestClose_full', 'MuscimolPostRestClose_cortex',...
    % 'SalinePreDecision_full', 'SalinePostDecision_cortex', 'SalinePreInfo_full', 'SalinePostInfo_cortex', ...
    % 'SalinePreInfoAnti_full', 'SalinePostInfoAnti_cortex','SalinePreInfoResp_full', 'SalinePostInfoResp_cortex',...
    % 'SalinePreRestOpen_full', 'SalinePostRestOpen_cortex','SalinePreRestClose_full', 'SalinePostRestClose_cortex',...
    % 'SimRecPreDecision_full', 'SimRecPreInfo_full',...
    % 'SimRecPreInfoAnti_full', 'SimRecPreInfoResp_full',...
    % 'SimRecPreRestOpen_full', 'SimRecPreRestClose_full',...
    % };
% task_names = {...
%     'MuscimolPreDecision_full', 'MuscimolPostDecision_cortex', 'MuscimolPreInfo_full', 'MuscimolPostInfo_cortex', ...
%     };
force_rebuild = true;
force_retrain = true;
debug = true;
total_training = 0;
skipped = 0;
failed= 0;
success = 0;

failed_list = {};

% task_names = {'MuscimolPreDecision_full', 'SalinePreDecision_full', 'SimRecPreDecision_full'};
% task_names = {'MuscimolPostDecision_full', 'SalinePostDecision_full'};
% trial_names = {'100B', '50BI', '50BN', '100S', '0'};
total_time = 0;
for task_idx=1:length(task_names)
    % for cuetype=1:5
    % compare if is Muscimol sessions
    if strcmp(task_names{task_idx}(1:8), 'Muscimol')
        session_idxs = [6,7,8,9,10,1,4,5,2,3];
    else
        session_idxs = [1,4,5,2,3];
    end
    % session_idxs = [1];

    for session_idx = session_idxs
        for trial_idx = 1:1
            try
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
                
                % task
                % session = 1;
                % kernel_name = 'expDecay10';
                % kernel_name = 'expMulti200';
                kernel_name = 'DeltaPure';
                
                % reg.l1=5;
                % reg.l2=0;
                % reg.name='L1=5';    
                
                reg.l1=0;
                reg.l2=2;
                reg.name='L2=2';
                
                % reg.l1=2;
                % reg.l2=0;
                % reg.name='L1=2';
                
                % reg.l1=0;
                % reg.l2=0;
                % reg.name='NoReg';
                
                shuffle_size=2;
                max_epoch=2500;
                
                
                %% generate shuffled raster
                fprintf("Shuffle rasters\n");
                tic;
                for shuffle_seed=1:shuffle_size
                    % skip if already exists
                    target_path = ['../GLM_data/', dataset_name, '/raster_', ...
                        dataset_name, '_', int2str(session_idx),  '_', int2str(shuffle_seed), '.mat'];
                    if isfile(target_path) && ~force_rebuild
                        fprintf("Skip %d. \n", shuffle_seed);
                        continue;
                    end

                    skip_flag = false;
                    shuffle_across_trial=(shuffle_seed<2);
                    shuffle(dataset_name, session_idx, shuffle_seed, shuffle_across_trial);
                end
                toc;

                %% convolve predj and combine trials
                fprintf("Convolution\n");
                tic;
                for shuffle_seed=0:shuffle_size % seed=0: original data (no shuffle)
                    % skip if already exists
                    target_path = ['../GLM_data/', dataset_name, '/GLMdata_', dataset_name,...
                        '_', int2str(session_idx), '_', kernel_name,  '_', int2str(shuffle_seed), '.mat'];
                    if isfile(target_path) && ~force_rebuild
                        fprintf("Skip %d. \n", shuffle_seed);
                        continue;
                    end

                    skip_flag = false;
                    convolution(dataset_name, session_idx, kernel_name, shuffle_seed);
                end
                toc;
                %% GLM inference
                for shuffle_seed=0:shuffle_size
                    fprintf("Training %d\n", shuffle_seed);

                    % skip if already exists
                    foldername = ['../GLM_model/', dataset_name];
                    target_path = [foldername, '/GLM_', dataset_name, '_', ...
                    int2str(session_idx), '_', kernel_name, '_', int2str(shuffle_seed), '_', ...
                    reg.name, '_', int2str(max_epoch)];
                    if isfile([target_path, '.mat']) && ~force_retrain
                        fprintf("Skip. \n");
                        continue;
                    end

                    skip_flag = false;
                    tic;
                    GLM_multi_kernel_err(dataset_name, session_idx, kernel_name, shuffle_seed, max_epoch, reg, 1, 5e-3);
                    toc;
                end

                %% plot
                % plot_GLM(dataset_name, session_idx, kernel_name, max_epoch, reg, shuffle_size);
                % type_file = ['../GLM_data/', dataset_name, '/celltype_', dataset_name, '_', int2str(session_idx), ...
                % '.mat'];
                % load(type_file, "cell_type");

                fprintf("Plotting\n");
                tic;
                channel_file = ['../GLM_data/', dataset_name, '/raster_', dataset_name, '_', int2str(session_idx), ...
                '_0.mat'];
                load(channel_file, "channel");
                plot_GLM_sorted(dataset_name, session_idx, kernel_name, max_epoch, reg, shuffle_size, "idx", channel);
                toc;

                %% plot gen
                % plot_generated(dataset_name)
                session_time = toc(tick_session);
                fprintf("Session time: %f\n", session_time);
                total_time = total_time + session_time;
                fprintf("Total time: %f\n", total_time);

                if skip_flag
                    skipped = skipped + 1;
                else
                    success = success + 1;
                end

            catch ME
                fprintf("Failed: %s\n", ME.message);
                failed = failed + 1;
                failed_list{end+1} = {dataset_name, int2str(session_idx), ME.message};
                if debug
                    throw(ME);
                end
            end
        end
    end
    % plot_rest(session_idx, kernel_name, max_epoch, reg, shuffle_size);
end

fprintf("Total: %d, Success: %d, Skipped: %d, Failed: %d\n", total_training, success, skipped, failed);
% save failed_list
save('../GLM_data/failed_list.mat', 'failed_list');
if failed>0
    for i=1:length(failed_list)
        fprintf("Failed: %s, %s, %s\n", failed_list{i}{1}, failed_list{i}{2}, failed_list{i}{3});
    end
end
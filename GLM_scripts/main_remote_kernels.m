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
%     'SalinePreDecision_cortex', 'SalinePostDecision_cortex', 'SalinePreInfo_cortex', 'SalinePostInfo_cortex', ...
%     'SalinePreInfoAnti_cortex', 'SalinePostInfoAnti_cortex','SalinePreInfoResp_cortex', 'SalinePostInfoResp_cortex',...
%     };
task_names = {'MuscimolPreDecision_full', 'SalinePreDecision_full', 'SimRecPreDecision_full'};
% trial_names = {'100B', '50BI', '50BN', '100S', '0'};
for task_idx=1:3
    % for cuetype=1:5
    if task_idx == 1
        session_idxs = [6,7,8,9,10,1,4,5,2,3];
    else
        session_idxs = [1,4,5,2,3];
    end

    for session_idx = session_idxs
        kernels = {'exp5Gauss5C20', 'exp5Gauss5C30', 'exp5Gauss5C40'};
        for kernel_idx = 1:3
            % if task_idx==1 && session_idx<10
            %     continue;
            % end
            % session_idx = session*100+cuetype;
            % fprintf("Main: session%d, cue%d\n", session, cuetype);

            % dataset_name = [task_names{task_idx}, '_', trial_names{trial_idx}];
            tick_session = tic;
            dataset_name = task_names{task_idx};
            fprintf("Main: %s, session%d\n", dataset_name, session_idx);

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
            % kernel_name = 'expGauss60';
            kernel_name = kernels{kernel_idx};
            
            % reg.l1=5;
            % reg.l2=0;
            % reg.name='L1=5';    
            
            reg.l1=0;
            reg.l2=2;
            reg.name='L2=2';
            
            % reg.l1=0;
            % reg.l2=0;
            % reg.name='NoReg';
            
            shuffle_size=2;
            max_epoch=2500;
            
            
            %% generate shuffled raster
            fprintf("Shuffle rasters\n");
            tic;
            for seed=1:shuffle_size
                shuffle_across_trial=(seed<2);
                shuffle(dataset_name, session_idx, seed, shuffle_across_trial);
            end
            toc;

            %% convolve predj and combine trials
            fprintf("Convolution\n");
            tic;
            for seed=0:shuffle_size % seed=0: original data (no shuffle)
                convolution(dataset_name, session_idx, kernel_name, seed);
            end
            toc;
            %% GLM inference
            for seed=0:shuffle_size
                tic;
                fprintf("Training %d\n", seed);
                GLM_multi_kernel(dataset_name, session_idx, kernel_name, seed, max_epoch, reg,1)
                toc;
            end
            % plot
            % plot_GLM(dataset_name, session_idx, kernel_name, max_epoch, reg, shuffle_size);
            % type_file = ['../GLM_data/', dataset_name, '/celltype_', dataset_name, '_', int2str(session_idx), ...
            % '.mat'];
            % load(type_file, "cell_type");
            fprintf("Plotting %d\n", seed);
            tic;
            channel_file = ['../GLM_data/', dataset_name, '/raster_', dataset_name, '_', int2str(session_idx), ...
            '_0.mat'];
            load(channel_file, "channel");
            plot_GLM_sorted(dataset_name, session_idx, kernel_name, max_epoch, reg, shuffle_size, "idx", channel);
            toc;
            %% plot gen
            % plot_generated(dataset_name)
            fprintf("Session time: %f\n", toc(tick_session));
        end
    end
    % plot_rest(session_idx, kernel_name, max_epoch, reg, shuffle_size);
end

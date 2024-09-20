%% GLM raster fitting
% Input: raw raster file "../GLM_data/[dataset_name]/raster_[dataset_name]_[session]_0.mat"

% Required variables in the input file:
% "n_trial": integer, trial number.
% "trial_len": int(1, n_trial), time bin number of each trial.
% "rasters": cell(1, n_trial), each element is a trial.
% Each raster is a (N, trial_len(i)) binary matrix, N is number of neurons.

clear
task_names = {'MuscimolPre_cortex', 'MuscimolPost_cortex', 'MuscimolPre_full'};
% task_names = {'MuscimolPostVLPFC'};
trial_names = {'100B', '50BI', '50BN', '100S', '0'};
for task_idx=3:3
    % for cuetype=1:5
    for session_idx=[1,2,3,4,5]
        for trial_idx=1:5
            % if task_idx==1 && session_idx<10
            %     continue;
            % end
            % session_idx = session*100+cuetype;
            % fprintf("Main: session%d, cue%d\n", session, cuetype);

            dataset_name = [task_names{task_idx}, '_', trial_names{trial_idx}];

            % dataset_name = task_names{task_idx};
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
            kernel_name = 'expDecay10';
            % kernel_name = 'steps25';
            
            % reg.l1=1;
            % reg.l2=0;
            % reg.name='L1=1';    
            
            reg.l1=0;
            reg.l2=1;
            reg.name='L2=1';
            
            % reg.l1=0;
            % reg.l2=0;
            % reg.name='NoReg';
            
            shuffle_size=6;
            max_epoch=2500;
            
            
            %% generate shuffled raster
            fprintf("Shuffle rasters\n");
            for seed=1:shuffle_size
                shuffle_across_trial=(seed<4);
                shuffle(dataset_name, session_idx, seed, shuffle_across_trial);
            end

            %% convolve predj and combine trials
            fprintf("Convolution\n");
            for seed=0:shuffle_size % seed=0: original data (no shuffle)
                convolution(dataset_name, session_idx, kernel_name, seed);
            end
            %% GLM inference
            for seed=0:shuffle_size
                fprintf("Training %d\n", seed);
                GLM_multi_kernel(dataset_name, session_idx, kernel_name, seed, max_epoch, reg,1)
            end
            %% plot
            plot_GLM(dataset_name, session_idx, kernel_name, max_epoch, reg, shuffle_size);
            % type_file = ['../GLM_data/', dataset_name, '/celltype_', dataset_name, '_', int2str(session_idx), ...
            % '.mat'];
            % load(type_file, "cell_type");
            channel_file = ['../GLM_data/', dataset_name, '/raster_', dataset_name, '_', int2str(session_idx), ...
            '_0.mat'];
            load(channel_file, "channel");
            plot_GLM_sorted(dataset_name, session_idx, kernel_name, max_epoch, reg, shuffle_size, "idx", channel);
            %% plot gen
            % plot_generated(dataset_name)
        end
    end
    % plot_rest(session_idx, kernel_name, max_epoch, reg, shuffle_size);
end

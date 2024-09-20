%% GLM raster fitting
% Input: raw raster file "../GLM_data/[dataset_name]/raster_[dataset_name]_[session]_0.mat"

% Required variables in the input file:
% "n_trial": integer, trial number.
% "trial_len": int(1, n_trial), time bin number of each trial.
% "rasters": cell(1, n_trial), each element is a trial.
% Each raster is a (N, trial_len(i)) binary matrix, N is number of neurons.

clear
task_names = {'eyeClose', 'eyeOpen'};
for session=1:5
    % for cuetype=1:5
    for sub_idx=1:2
        % session_idx = session*100+cuetype;
        % fprintf("Main: session%d, cue%d\n", session, cuetype);
        fprintf("Main: session%d, %s\n", session, task_names{sub_idx});
        session_idx = session;
    %% parameters
        % dataset_name = 'generated';
        % % session = 0;
        % kernel_name = 'expDecay10';
        
        % dataset_name = 'taskCue';
        % % session = 1;
        % kernel_name = 'expDecay10';
        
        % eyeclose/eyeopen
        dataset_name = task_names{sub_idx};
        % session = 1;
        kernel_name = 'expDecay10';
        
        reg.l1=5;
        reg.l2=0;
        reg.name='L1=5';    
        
        % reg.l1=0;
        % reg.l2=1;
        % reg.name='L2=1';
        
        % reg.l1=0;
        % reg.l2=0;
        % reg.name='NoReg';
        
        shuffle_size=8;
        max_epoch=2500;
        
        
        %% generate shuffled raster
        % fprintf("Shuffle rasters\n");
        % for seed=1:shuffle_size
        %     shuffle(dataset_name, session_idx, seed);
        % end

        % %% convolve predj and combine trials
        % fprintf("Convolution\n");
        % for seed=0:shuffle_size % seed=0: original data (no shuffle)
        %     convolution(dataset_name, session_idx, kernel_name, seed);
        % end
        % %% GLM inference
        % for seed=0:shuffle_size
        %     fprintf("Training %d\n", seed);
        %     GLM_multi_kernel(dataset_name, session_idx, kernel_name, seed, max_epoch, reg,1)
        % end

        %% plot
        % plot_GLM(dataset_name, session_idx, kernel_name, max_epoch, reg, shuffle_size);
        % plot_GLM_sorted(dataset_name, session_idx, kernel_name, max_epoch, reg, shuffle_size, "max");
        %% plot gen
        % plot_generated(dataset_name)
    end
    plot_rest(session_idx, kernel_name, max_epoch, reg, shuffle_size);
end

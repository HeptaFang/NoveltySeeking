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
all_diffs = cell(3, 3);
all_Js = cell(3, 3, 2);
all_Js_abs = cell(3, 3, 2);
kernel_diffs = cell(3, 3, 3);
kernel_Js = cell(3, 3, 2, 3);
kernel_Js_abs = cell(3, 3, 2, 3);

task_names = {'MuscimolPreDecision_full', 'SalinePreDecision_full', 'SimRecPreDecision_full'};
% trial_names = {'100B', '50BI', '50BN', '100S', '0'};
for task_idx=1:3
    % for cuetype=1:5
    if task_idx < 2
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
            
            
            % %% generate shuffled raster
            % fprintf("Shuffle rasters\n");
            % tic;
            % for seed=1:shuffle_size
            %     shuffle_across_trial=(seed<3);
            %     shuffle(dataset_name, session_idx, seed, shuffle_across_trial);
            % end
            % toc;

            % %% convolve predj and combine trials
            % fprintf("Convolution\n");
            % tic;
            % for seed=0:shuffle_size % seed=0: original data (no shuffle)
            %     convolution_gpu(dataset_name, session_idx, kernel_name, seed);
            % end
            % toc;
            % %% GLM inference
            % for seed=0:shuffle_size
            %     tic;
            %     fprintf("Training %d\n", seed);
            %     GLM_multi_kernel(dataset_name, session_idx, kernel_name, seed, max_epoch, reg,1)
            %     toc;
            % end
            % % plot
            % % plot_GLM(dataset_name, session_idx, kernel_name, max_epoch, reg, shuffle_size);
            % % type_file = ['../GLM_data/', dataset_name, '/celltype_', dataset_name, '_', int2str(session_idx), ...
            % % '.mat'];
            % % load(type_file, "cell_type");
            % fprintf("Plotting %d\n", seed);

            tic;
            channel_file = ['../GLM_data/', dataset_name, '/raster_', dataset_name, '_', int2str(session_idx), ...
            '_0.mat'];
            load(channel_file, "channel");
            [session_diff, session_J, session_J_abs] = kernel_weight_comparason(dataset_name, session_idx, kernel_name, max_epoch, reg, shuffle_size, "idx", channel);
            for i = 1:3
                for j = 1:3
                    all_diffs{i, j} = [all_diffs{i, j}; session_diff{i, j}];
                    all_Js{i, j, 1} = [all_Js{i, j, 1}; session_J{i, j, 1}];
                    all_Js{i, j, 2} = [all_Js{i, j, 2}; session_J{i, j, 2}];
                    all_Js_abs{i, j, 1} = [all_Js_abs{i, j, 1}; session_J_abs{i, j, 1}];
                    all_Js_abs{i, j, 2} = [all_Js_abs{i, j, 2}; session_J_abs{i, j, 2}];
                    kernel_diffs{i, j, kernel_idx} = [kernel_diffs{i, j, kernel_idx}; session_diff{i, j}];
                    kernel_Js{i, j, 1, kernel_idx} = [kernel_Js{i, j, 1, kernel_idx}; session_J{i, j, 1}];
                    kernel_Js{i, j, 2, kernel_idx} = [kernel_Js{i, j, 2, kernel_idx}; session_J{i, j, 2}];
                    kernel_Js_abs{i, j, 1, kernel_idx} = [kernel_Js_abs{i, j, 1, kernel_idx}; session_J_abs{i, j, 1}];
                    kernel_Js_abs{i, j, 2, kernel_idx} = [kernel_Js_abs{i, j, 2, kernel_idx}; session_J_abs{i, j, 2}];
                end
            end
            toc;

            %% plot gen
            % plot_generated(dataset_name)
        end
    end
    % plot_rest(session_idx, kernel_name, max_epoch, reg, shuffle_size);
end

% load("../GLM_data/all_diffs.mat", "all_diffs");
% load("../GLM_data/kernel_diffs.mat", "kernel_diffs");
save("../GLM_data/all_diffs.mat", "all_diffs");
save("../GLM_data/all_Js.mat", "all_Js");
save("../GLM_data/all_Js_abs.mat", "all_Js_abs");
save("../GLM_data/kernel_diffs.mat", "kernel_diffs");
save("../GLM_data/kernel_Js.mat", "kernel_Js");
save("../GLM_data/kernel_Js_abs.mat", "kernel_Js_abs");

fprintf("Plotting all sessions\n");
tic;
kernel_weight_comparason_allsession(all_diffs, 'all');
for kernel_idx = 1:3
    kernel_weight_comparason_allsession(kernel_diffs(:, :, kernel_idx), kernels{kernel_idx});
end
toc;
function calc_pred(dataset_name, session_idx, kernel_name, epoch, reg)
% Generate predicted p from predj and model parameters
% Integrated into inference.

%% Load predj
file_path = ['../GLM_data/', session_stage_full,'/GLMdata_',session_stage_full, '_', ...
                app.session_idx,'_', app.kernel, '_0.mat'];



%% Load model
file_path = ['../GLM_model/', session_stage_full,...
                    '/GLM_', session_stage_full, '_', app.session_idx, '_',...
                    app.kernel, '_0_', app.reg, '_', app.epoch, '.mat'];




end
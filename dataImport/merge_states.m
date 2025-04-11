% Merge states: task, eyes open, and eyes closed

session_types = {'Muscimol', 'Saline'};
states = {'Task', 'RestOpen', 'RestClose'};
areas = {'ACC', 'Thalamus', 'VLPFC'};
n_states = length(states);
n_areas = length(areas);

kernel = 'DeltaPure';
align = 'AlignLast';

for session_type_idx = 1:2
    session_type = session_types{session_type_idx};
    if strcmp(session_type, 'Muscimol')
        n_session = 10;
    else
        n_session = 5;
    end

    for session_idx = 1:n_session
        for prepost = 1:3

            %% merge GLMdata
            B_total = 0;
            predjs_PS_total = NaN;
            predjs_conn_total = NaN;
            raster_total = NaN;

            for state_idx = 1:n_states
                state = states{state_idx};
                
                switch prepost
                    case 1
                        session_stage = [session_type, 'Pre', state, '_full'];
                    case 2
                        session_stage = [session_type, 'Post', state, '_cortex'];
                    case 3
                        session_stage = [session_type, 'Pre', state, '_cortex'];
                end
                fprintf('Loading %s %s %s\n', session_type, state, int2str(session_idx));
                
                % Load data
                file_path = ['../GLM_data/', session_stage, '_', align, '/GLMdata_', session_stage,...
                 '_', align, '_', int2str(session_idx), '_', kernel, '_0.mat'];
                load(file_path, 'N', 'B', 'PS_kernels', 'conn_kernels', 'kernel_len',...
                 'n_PS_kernel', 'n_conn_kernel', 'predjs_PS', 'predjs_conn', 'raster');
                B_total = B_total + B;

                if all(isnan(predjs_PS_total))
                    predjs_PS_total = predjs_PS;
                    predjs_conn_total = predjs_conn;
                    raster_total = raster;
                else
                    predjs_PS_total = cat(2, predjs_PS_total, predjs_PS);
                    predjs_conn_total = cat(2, predjs_conn_total, predjs_conn);
                    raster_total = cat(2, raster_total, raster);
                end
            end

            % Save data
            switch prepost
                case 1
                    session_stage_new = [session_type, 'Pre', 'All', '_full'];
                case 2
                    session_stage_new = [session_type, 'Post', 'All', '_cortex'];
                case 3
                    session_stage_new = [session_type, 'Pre', 'All', '_cortex'];
            end

            file_folder = ['../GLM_data/', session_stage_new, '_', align, '/'];
            if ~exist(file_folder, 'dir')
                mkdir(file_folder);
            end
            file_path = [file_folder, 'GLMdata_', session_stage_new, '_', align, '_', int2str(session_idx), '_', kernel, '_0.mat'];
            B = B_total;
            predjs_conn = predjs_conn_total;
            predjs_PS = predjs_PS_total;
            raster = raster_total;
            save(file_path, 'N', 'B', 'PS_kernels', 'conn_kernels', 'kernel_len',...
             'n_PS_kernel', 'n_conn_kernel', 'predjs_PS', 'predjs_conn', 'raster');
            fprintf('Saved %s %s\n', session_type, int2str(session_idx));

            %% copy borders and rasters
            file_path = ['../GLM_data/', session_stage, '_', align, '/raster_', session_stage, '_', align, '_', int2str(session_idx), '_0.mat'];
            % load(file_path, 'n_trial', 'cell_id', 'cell_area', 'channel', 'session_name_full', 'trial_len', 'rasters');
            file_path_new = ['../GLM_data/', session_stage_new, '_', align, '/raster_', session_stage_new, '_', align, '_', int2str(session_idx), '_0.mat'];
            copyfile(file_path, file_path_new);

            file_path = ['../GLM_data/',  session_stage, '_', align,  '/borders_', session_stage, '_', align, '_', int2str(session_idx), '.mat'];
            % load(file_path, 'borders');
            file_path_new = ['../GLM_data/', session_stage_new, '_', align, '/borders_', session_stage_new, '_', align, '_', int2str(session_idx), '.mat'];
            copyfile(file_path, file_path_new);
        end
    end
end
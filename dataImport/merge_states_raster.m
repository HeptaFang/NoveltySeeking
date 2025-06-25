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

            %% merge rasters
            n_trial_total = 0;
            rasters_all = cell(1, 0);
            firing_rates_all = cell(1, 0);
            trial_len_all = [];

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
                file_path = ['../GLM_data/', session_stage, '_', align, '/raster_',... 
                session_stage, '_', align, '_', int2str(session_idx), '_0.mat'];
                load(file_path, 'cell_area', 'cell_id', 'channel', 'firing_rates',...
                    'n_trial', 'rasters', 'session_name_full', 'trial_len');
                n_trial_total = n_trial_total + n_trial;
                rasters_all = [rasters_all, rasters];
                firing_rates_all = [firing_rates_all, firing_rates];
                trial_len_all = [trial_len_all, trial_len];
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
            file_path_new = ['../GLM_data/', session_stage_new, '_', align, '/raster_',...
                session_stage_new, '_', align, '_', int2str(session_idx), '_0.mat'];
            
            n_trial = n_trial_total;
            rasters = rasters_all;
            firing_rates = firing_rates_all;
            trial_len = trial_len_all;
            save(file_path_new, 'cell_area', 'cell_id', 'channel', 'firing_rates',...
             'n_trial', 'rasters', 'session_name_full', 'trial_len');
            fprintf('Saved %s %s\n', session_type, int2str(session_idx));

            %% copy borders
            file_path = ['../GLM_data/',  session_stage, '_', align,  '/borders_', session_stage, '_', align, '_', int2str(session_idx), '.mat'];
            file_path_new = ['../GLM_data/', session_stage_new, '_', align, '/borders_', session_stage_new, '_', align, '_', int2str(session_idx), '.mat'];
            copyfile(file_path, file_path_new);
        end
    end
end
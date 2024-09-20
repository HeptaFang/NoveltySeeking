
% load data
tasknames = {'MuscimolPre_full','MuscimolPre_cortex','MuscimolPost_cortex'};
for task_idx = 3:3
    taskname = tasknames{task_idx};
    for session = 1:10
        folder = ['../GLM_data/', taskname];
        filename = [folder, '/raster_', taskname, '_', int2str(session), '_0.mat'];
        full = load(filename);
        border_file = [folder, '/borders_', taskname, '_', int2str(session), '.mat'];
        load(border_file, "borders");
        channel_file = ['../GLM_data/', taskname, '/raster_', taskname, '_', int2str(session), ...
        '_0.mat'];
        load(channel_file, "channel");

        cuetypes = [6710, 6720, 6730, 6740, 6750];
        cue_names = {'100B', '50BI', '50BN', '100S', '0'};
        for i=1:5
            mask = (full.cuetype == cuetypes(i));
            
            N = full.N;
            cell_area = full.cell_area;
            cell_id = full.cell_id;
            cuetype = cuetypes(i);
            firing_rates = full.firing_rates(mask);
            n_trial = sum(mask);
            rasters = full.rasters(mask);
            session_name_full = full.session_name_full;
            spikes = full.spikes(:, mask);
            trial_len = full.trial_len(mask);

            taskname_new = [taskname, '_', cue_names{i}];
            check_path(['../GLM_data/', taskname_new]);
            savefile = ['../GLM_data/', taskname_new, '/raster_',...
                taskname_new, '_', int2str(session), '_0.mat'];
            save(savefile, "N", "cell_id", "cuetype", "trial_len","spikes",...
                "session_name_full","rasters","n_trial","firing_rates", "cell_area", "channel");
            border_file = ['../GLM_data/', taskname_new,'/borders_', taskname_new, '_', int2str(session), '.mat'];
            save(border_file, "borders");
        end
    end
end
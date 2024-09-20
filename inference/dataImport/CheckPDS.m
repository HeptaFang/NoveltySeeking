% load cell-based PDS file, 
% extract trial-based spike timings and rasters.
% 

% 1ms time bin, max time = 3.1s
dt=0.001;
B=3100;
max_t=B*dt;
edges = 0:dt:max_t;

unique_sessions_all = ...
    {{'10272023', '11012023', '11102023', '11172023', '12012023',...
    '12082023', '12152023', '12292023', '01052024', '01122024'},...
    {'01302024', '02022024', '02092024', '02162024', '02292024'}};

% load data
controls = {'Muscimol', 'Saline'};
areas = {'ACC', 'Thalamus', 'VLPFC'};
for control_idx = 1:2
    control = controls{control_idx};
    unique_sessions = unique_sessions_all{control_idx};
    session_num = length(unique_sessions);

    for area_idx = 1:3
        area = areas{area_idx};
        folder_name = ['../MYInfoPavChoice/', control, '/', area];
        filenames = {dir(folder_name).name}.';
        splited = cellfun(@(name) split(name, ["-", "_"]), filenames, 'UniformOutput', false);
        splited = splited(3:end);

        session_names = cellfun(@(x) x{2}, splited, 'UniformOutput', false);
        session_types = cellfun(@(x) x{3}, splited, 'UniformOutput', false);
        cell_ids = cellfun(@(x) x{5}, splited, 'UniformOutput', false);


        for session_idx = 1:session_num
            % metadata for sessions, need to adjust for each dataset.
            session_name = unique_sessions{session_idx};
            session_type = '008';
            if strcmp(area, 'Thalamus')
                session_type = '001';
                if strcmp(session_name, '11012023')
                    session_type = '004';
                end
            end
            
            session_file_idx = strcmp(session_names, session_name)...
                & strcmp(session_types, session_type);
            N = sum(session_file_idx);
            cell_id = cell_ids(session_file_idx).';

            % task trials
            subsession_types = {'Pre', 'Post'};
            taskIDs = [1.1, 1.2];
            for subsession_idx = 1:2
                subsession = subsession_types{subsession_idx};
                taskID = taskIDs(subsession_idx);

                trial_num = NaN;
                spikes = NaN;
                rasters = NaN;
                firing_rates = NaN;
                cuetype = NaN;
                
                % construct output file
                session = struct();
                session_name_full = [session_name, '_', subsession, '_', area];
                session.N = N;

                fprintf("Loading: %s, %s, %s, session%d, N=%d...\n", control, area, subsession, session_idx, N);

                for i = 1:N
                    % fprintf("Loading: %s, %s, %s, session%d, #%d...\n", control, area, subsession, session_idx, i);
                    filename = [folder_name, '/LemmyKim-', session_name, '-', session_type, ...
                        '_MYInfoPavChoice_', cell_id{i}, '_PDS.mat'];
                    load(filename, "PDS");
                    fprintf("%d.\tcell ID: %s,\ttrial=%d\n", i, cell_id{i}, size(PDS.trialnumber, 2));
                    
                end

            end
            % todo: load resting state

        end
    end
end
% load cell-based PDS file, 
% extract trial-based spike timings and rasters.
% 

% 1ms time bin, max time = 3.1s
dt=0.001;
% B=3100;
% max_t=B*dt;
% edges = 0:dt:max_t;

unique_sessions_all = ...
    {{'10272023', '11012023', '11102023', '11172023', '12012023',...
    '12082023', '12152023', '12292023', '01052024', '01122024'},...
    {'01302024', '02022024', '02092024', '02162024', '02292024'},...
    {'08112023', '08142023', '08152023', '08162023', '08172023'}};

% load data
controls = {'Muscimol', 'Saline', 'SimRec'};
areas = {'ACC', 'Thalamus', 'VLPFC'};
for control_idx = 1:3
    control = controls{control_idx};
    unique_sessions = unique_sessions_all{control_idx};
    session_num = length(unique_sessions);

    randomA_start = NaN;
    % randomA_end = NaN;
    randomB_start = NaN;
    % randomB_end = NaN;
    randomShort_start = NaN;
    % randomShort_end = NaN;
    randomLong_start = NaN;
    % randomLong_end = NaN;

    for area_idx = 1:3
        area = areas{area_idx};
        folder_name = ['../', control, '/', area];
        filenames = {dir(folder_name).name}.'; % format: LemmyKim-#date-00#-MYInfoPavChoice_cl##_PDS.mat
        splited = cellfun(@(name) split(name, ["-", "_"]), filenames, 'UniformOutput', false);
        splited = splited(3:end);

        session_names = cellfun(@(x) x{2}, splited, 'UniformOutput', false);
        session_types = cellfun(@(x) x{3}, splited, 'UniformOutput', false);
        if strcmp(control, 'SimRec')
            cell_ids = cellfun(@(x) x{6}, splited, 'UniformOutput', false);
        else
            cell_ids = cellfun(@(x) x{5}, splited, 'UniformOutput', false);
        end

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
            if strcmp(control, 'SimRec')
                session_type = '001';
                if strcmp(session_name, '08112023')
                    session_type = '003';
                end
            end
            
            session_file_idx = strcmp(session_names, session_name)...
                & strcmp(session_types, session_type);
            N = sum(session_file_idx);
            cell_id = cell_ids(session_file_idx).';

            % trialphases = {'Decision', 'InfoAnti', 'InfoResp', 'Info', 'RandomA', 'RandomB'};
            % trialphases = {'Offer1', 'Offer2', 'Decision', 'InfoAnti', 'InfoResp', 'Reward', 'RandomA', 'RandomB'};
            trialphases = {'Task', 'RandomShort', 'RandomLong'};
            for trialphase_idx = 1:3
                trialphase = trialphases{trialphase_idx};

                % task trials
                subsession_types = {'Pre', 'Post'};
                
                % simrec and other two have different taskIDs
                if strcmp(control, 'SimRec')
                    taskIDs = [1, 2];
                else
                    taskIDs = [1.1, 1.2];
                end
                for subsession_idx = 1:2
                    subsession = subsession_types{subsession_idx};
                    taskID = taskIDs(subsession_idx);

                    % skip simrec post sessions
                    if strcmp(control, 'SimRec') && strcmp(subsession, 'Post')
                        continue;
                    end
                    % % skip muscimol thalamus post sessions
                    % if strcmp(area, 'Thalamus') && strcmp(subsession, 'Post') && strcmp(control, 'Muscimol')
                    %     continue;
                    % end
                    % skip all thalamus post sessions
                    if strcmp(area, 'Thalamus') && strcmp(subsession, 'Post')
                        continue;
                    end

                    % Loading a single file starts here

                    trial_num = NaN;
                    spikes = NaN;
                    rasters = NaN;
                    firing_rates = NaN;
                    cuetype = NaN;
                    channel = NaN;
                    trial_len = NaN;
                    
                    % construct output file
                    session_name_full = [session_name, '_', subsession, trialphase, '_', area];

                    fprintf("Loading: %s, %s, %s, %s, session%d, N=%d...\n", control, area, subsession, trialphase, session_idx, N);
                    max_duration = 0;
                    min_duration = 9999;
                    max_trial_len = 0;
                    min_trial_len = 999999;
                    for i = 1:N
                        % fprintf("Loading: %s, %s, %s, session%d, #%d...\n", control, area, subsession, session_idx, i);
                        if strcmp(control, 'SimRec')
                            filename = [folder_name, '/LemmyKim-', session_name, '-', session_type, ...
                                '_MYInfoPavChoice_NovelReward_', cell_id{i}, '_PDS.mat'];
                        else
                            filename = [folder_name, '/LemmyKim-', session_name, '-', session_type, ...
                                '_MYInfoPavChoice_', cell_id{i}, '_PDS.mat'];
                        end
                        load(filename, "PDS");
                        
                        % NaN check for thalamus
                        if any(isnan(PDS.taskID))
                            % fprintf("NaN in: %s, %s, %s, session%d, #%d...\n", control, area, subsession, session_idx, i);
                            PDS.taskID(isnan(PDS.taskID)) = 1.1;
                        end
                        
                        % ---choose trials
                        % -choice and pav trials
                        % selected_trial = (PDS.taskID==taskID)&(PDS.goodtrial);
                        % -pav trials only
                        % selected_trial = (PDS.taskID==taskID)&(PDS.goodtrial)&isnan(PDS.timeoffer1cho)&isnan(PDS.timeoffer2cho);
                        % -choice trials only
                        selected_trial = (PDS.taskID==taskID)&(PDS.goodtrial)&(~isnan(PDS.timeoffer1cho)|~isnan(PDS.timeoffer2cho));

                        selected_spikes = PDS.sptimes(selected_trial);

                        % ---choose phases in trial
                        % calc decision time
                        decision_time = [PDS.timeoffer1cho;PDS.timeoffer2cho];
                        decision_time(isnan(decision_time)) = 0;
                        decision_time = sum(decision_time, 1);

                        switch trialphase
                            case 'All'
                            % -all
                            selected_start = PDS.timetargeton(selected_trial);
                            selected_end = PDS.timeInfotargetoff(selected_trial);

                            case 'Offer1'
                            % -after offer1
                            selected_start = PDS.timetargeton(selected_trial);
                            selected_end = selected_start + 1;

                            case 'Offer2'
                            % -after offer2
                            selected_start = PDS.timetargeton1(selected_trial);
                            selected_end = selected_start + 1;

                            case 'Decision'
                            % -before decision
                            selected_end = decision_time(selected_trial);
                            selected_start = selected_end - 1;

                            case 'InfoAnti'
                            % -before info cue
                            selected_end = PDS.timeInfotargeton(selected_trial);
                            selected_start = selected_end - 1;

                            case 'InfoResp'
                            % -after info cue
                            selected_start = PDS.timeInfotargeton(selected_trial);
                            selected_end = selected_start + 1;

                            case 'Reward'
                            % -after reward
                            selected_start = PDS.timereward(selected_trial);
                            selected_end = selected_start + 1;

                            case 'RandomA'
                            % -random 1s in the trial
                            if ~isa(randomA_start, "cell")
                                randomA_start = cell(session_num, 2);
                            end
                            if isempty(randomA_start{session_idx, subsession_idx})
                                randomA_start{session_idx, subsession_idx} = rand(1, sum(selected_trial));
                            end
                            trial_range = PDS.trialendtime(selected_trial) - PDS.trialstarttime(selected_trial) - 1;
                            selected_start = randomA_start{session_idx, subsession_idx} .* trial_range;
                            selected_end = selected_start + 1;

                            case 'RandomB'
                            % -random 1s in the trial
                            if ~isa(randomB_start, "cell")
                                randomB_start = cell(session_num, 2);
                            end
                            if isempty(randomB_start{session_idx, subsession_idx})
                                randomB_start{session_idx, subsession_idx} = rand(1, sum(selected_trial));
                            end
                            trial_range = PDS.trialendtime(selected_trial) - PDS.trialstarttime(selected_trial) - 1;
                            selected_start = randomB_start{session_idx, subsession_idx} .* trial_range;
                            selected_end = selected_start + 1;

                            case 'RandomShort'
                            % -random 2s in the trial
                            if ~isa(randomShort_start, "cell")
                                randomShort_start = cell(session_num, 2);
                            end
                            if isempty(randomShort_start{session_idx, subsession_idx})
                                randomShort_start{session_idx, subsession_idx} = rand(1, sum(selected_trial));
                            end
                            trial_range = PDS.trialendtime(selected_trial) - PDS.trialstarttime(selected_trial) - 2.00005;
                            selected_start = randomShort_start{session_idx, subsession_idx} .* trial_range;
                            selected_end = selected_start + 2.00005;

                            case 'RandomLong'
                            % -random 5s in the trial
                            if ~isa(randomLong_start, "cell")
                                randomLong_start = cell(session_num, 2);
                            end
                            if isempty(randomLong_start{session_idx, subsession_idx})
                                randomLong_start{session_idx, subsession_idx} = rand(1, sum(selected_trial));
                            end
                            trial_range = PDS.trialendtime(selected_trial) - PDS.trialstarttime(selected_trial) - 5.00005;
                            selected_start = randomLong_start{session_idx, subsession_idx} .* trial_range;
                            selected_end = selected_start + 5.00005;

                            case 'Task'
                            % -full task
                            selected_start = PDS.trialstarttime(selected_trial);
                            selected_end = PDS.trialendtime(selected_trial);
                            
                        end

                        max_dur_cell = max(selected_end - selected_start);
                        min_dur_cell = min(selected_end - selected_start);
                        if max_dur_cell > max_duration
                            max_duration = max_dur_cell;
                        end
                        if min_dur_cell < min_duration
                            min_duration = min_dur_cell;
                        end

                        % if this is the first neuron, then initialize
                        if isnan(trial_num)
                            trial_num = sum(selected_trial);
                            spikes = cell(N, trial_num);
                            rasters = cell(1, trial_num);
                            firing_rates = cell(1, trial_num);
                            cuetype = PDS.Cuetype(selected_trial);
                            channel = zeros(1, N);
                            trial_len = ones(1, trial_num) * NaN;
                            for j=1:trial_num
                                rasters{j} = NaN;
                                firing_rates{j} = ones(N, 1) * NaN;
                            end
                        else
                            if trial_num~=sum(selected_trial)
                                error('trial num not match');
                            end
                        end

                        channel(1, i) = PDS.channel;
                        % spikes & rasters for each trial
                        for j=1:trial_num
                            spike_trial = selected_spikes{j};

                            % only keep spikes from cue to reward,
                            % align to cue.
                            spike_trial = spike_trial(spike_trial>selected_start(j) & spike_trial<selected_end(j));
                            spike_trial = spike_trial - selected_start(j);
                            
                            % auto fit start and end
                            trial_duration = selected_end(j) - selected_start(j);

                            % fixed trial duration
                            % trial_duration = 1;

                            trial_edges = 0:dt:trial_duration;
                            trial_B = length(trial_edges) - 1;
                            trial_len(j) = trial_B;

                            if trial_B > max_trial_len
                                max_trial_len = trial_B;
                            end
                            if trial_B < min_trial_len
                                min_trial_len = trial_B;
                            end

                            spikes{i, j} = spike_trial;
                            raster = histcounts(spike_trial, trial_edges);
                            raster(raster>1) = 1;
                            if isnan(rasters{j})
                                rasters{j} = zeros(N, trial_B);
                            end
                            rasters{j}(i, :) = raster;
                            firing_rates{j}(i) = mean(raster);
                        end
                    end

                    fprintf('max duration: %d, min duration: %d\n', max_duration, min_duration);
                    fprintf('max length: %d, min length: %d\n\n', ceil(max_duration/dt), ceil(min_duration/dt));
                    
                    % output
                    n_trial = trial_num;
                    % if isnan(n_trial)
                    %     trial_len = NaN;
                    % else
                    %     trial_len = ones(1, n_trial) * B;% fix this
                    % end

                    % save
                    save_folder = ['../GLM_data/', control, subsession, trialphase, '_', area];
                    check_path(save_folder);
                    save([save_folder, '/raster_', control, subsession, trialphase, '_', area,...
                        '_', int2str(session_idx), '_0.mat'],...
                        "rasters", "spikes", "firing_rates", "n_trial", "trial_len", ...
                        "session_name_full", "N", "cuetype", "cell_id", "channel");
                end
            % todo: load resting state
            end
        end
    end
end
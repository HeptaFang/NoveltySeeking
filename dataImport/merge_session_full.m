unique_sessions_all = ...
    {{'10272023', '11012023', '11102023', '11172023', '12012023',...
    '12082023', '12152023', '12292023', '01052024', '01122024'},...
    {'01302024', '02022024', '02092024', '02162024', '02292024'},...
    {'08112023', '08142023', '08152023', '08162023', '08172023'}};

% load data
controls = {'Muscimol', 'Saline', 'SimRec'};
% controls = {'Muscimol', 'Saline'};
areas = {'ACC', 'Thalamus', 'VLPFC'};
% areas = {'ACC', 'VLPFC'};
for control_idx = 1:3
    control = controls{control_idx};
    unique_sessions = unique_sessions_all{control_idx};
    session_num = length(unique_sessions);

    for session_idx = 1:session_num
        session_name = unique_sessions{session_idx};
        % subsession_types = {'Pre', 'Post'};
        subsession_types = {...
            'PreOffer1', 'PreOffer2','PreDecision', 'PreInfoAnti',...
            'PreInfoResp', 'PreReward', 'PreRandomA', 'PreRandomB',...
            'PostOffer1', 'PostOffer2','PostDecision', 'PostInfoAnti',...
            'PostInfoResp', 'PostReward', 'PostRandomA', 'PostRandomB',...
            };
        % subsession_types = {'PostDecision'};
        for subsession_idx = 1:length(subsession_types)
            subsession = subsession_types{subsession_idx};
            % skip muscimol and SimRec post sessions
            if (strcmp(control, 'Muscimol')||strcmp(control, 'SimRec')) && strcmp(subsession(1:4), 'Post')
                continue;
            end
            % skip SimRec RestClose and RestOpen
            if strcmp(control, 'SimRec') && (strcmp(subsession, 'PreRestClose') || strcmp(subsession, 'PreRestOpen'))
                continue;
            end

            fprintf('-------------------\n');
            fprintf('Merging: %s, %s, session%d...\n\n', control, subsession, session_idx);

            tic;
            N=0;
            cell_id = cell(1, 0);
            cuetype = zeros(1, 0);
            n_trial = NaN;
            cell_area = cell(1, 0);
            session_name_full = [session_name, '_', subsession];
            trial_len = NaN;
            % disp(session_name_full);

            borders = [];

            for area_idx = 1:length(areas)
                area = areas{area_idx};
                % load data from multiple areas
                
                area_file = ['../GLM_data/', control, subsession, '_', area, '/raster_', ...
                    control, subsession, '_', area, '_', int2str(session_idx), '_0.mat'];
                % check if data exists
                if ~isfile(area_file)
                    fprintf('Skipping %s: File not found.\n', area);
                    continue;
                end
                data = load(area_file);

                borders = [borders, N+data.N+0.5];

                if data.N==0
                    continue;
                end

                if isnan(n_trial)
                    n_trial = data.n_trial;
                    trial_len = data.trial_len;
                    cuetype = data.cuetype;
                    rasters = cell(1, n_trial);
                    firing_rates = cell(1, n_trial);
                    spikes = cell(0, n_trial);
                    channel = zeros(1, 0);
                    for i=1:n_trial
                        rasters{1, i} = zeros(0, trial_len(i));
                        firing_rates{1, i} = zeros(0, 1);
                    end
                else
                    if n_trial ~= data.n_trial
                        error('trial num not match! Area: %s, Control: %s, Subsession: %s', area, control, subsession);
                    end
                    if (any(cuetype~=data.cuetype & ~isnan(cuetype))) || any(trial_len ~= data.trial_len)
                        error('trial info not match! Area: %s, Control: %s, Subsession: %s', area, control, subsession);
                    end
                end
                N = N+data.N;
                cell_area = [cell_area, repmat({area}, 1, data.N)];
                cell_id = [cell_id, data.cell_id];
                spikes = [spikes;data.spikes];
                channel = [channel, data.channel];
                for i=1:n_trial
                    rasters{1, i} = [rasters{1, i}; data.rasters{1, i}];
                    firing_rates{1, i} = [firing_rates{1, i}; data.firing_rates{1, i}];
                end
            end
            % save data
            check_path(['../GLM_data/', control, subsession, '_full']);
            merged_file = ['../GLM_data/', control, subsession, '_full/raster_', ...
                    control, subsession, '_full_', int2str(session_idx), '_0.mat'];
            save(merged_file, 'N', "cell_id", "cuetype", "firing_rates", "n_trial",...
                "rasters", "spikes", "session_name_full", "trial_len", "cell_area", "channel");
            
            save(['../GLM_data/',control, subsession, '_full/borders_',...
                control, subsession, '_full_', ...
                int2str(session_idx),'.mat'], "borders");
            toc;
        end
    end
end
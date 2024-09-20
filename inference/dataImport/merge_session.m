unique_sessions_all = ...
    {{'10272023', '11012023', '11102023', '11172023', '12012023',...
    '12082023', '12152023', '12292023', '01052024', '01122024'},...
    {'01302024', '02022024', '02092024', '02162024', '02292024'}};

% load data
controls = {'Muscimol', 'Saline'};
% areas = {'ACC', 'Thalamus', 'VLPFC'};
areas = {'ACC', 'VLPFC'};
for control_idx = 1:2
    control = controls{control_idx};
    unique_sessions = unique_sessions_all{control_idx};
    session_num = length(unique_sessions);

    for session_idx = 1:session_num
        session_name = unique_sessions{session_idx};
        subsession_types = {'Pre', 'Post'};
        for subsession_idx = 1:2
            subsession = subsession_types{subsession_idx};

            N=0;
            cell_id = cell(1, 0);
            cuetype = zeros(1, 0);
            n_trial = NaN;
            cell_area = cell(1, 0);
            session_name_full = [session_name, '_', subsession];
            trial_len = NaN;

            for area_idx = 1:length(areas)
                area = areas{area_idx};
                % load data from multiple areas
                
                area_file = ['../GLM_data/', control, subsession, area, '/raster_', ...
                    control, subsession, area, '_', int2str(session_idx), '_0.mat'];
                data = load(area_file);

                if isnan(n_trial)
                    n_trial = data.n_trial;
                    trial_len = data.trial_len;
                    cuetype = data.cuetype;
                    rasters = cell(1, n_trial);
                    firing_rates = cell(1, n_trial);
                    spikes = cell(0, n_trial);
                    for i=1:n_trial
                        rasters{1, i} = zeros(0, trial_len);
                        firing_rates{1, i} = zeros(0, 1);
                    end
                else
                    if n_trial ~= data.n_trial || trial_len ~= data.trial_len
                        error('trial info not match!');
                    end
                    if any(cuetype~=data.cuetype)
                        error('trial type not match!');
                    end
                end
                N = N+data.N;
                cell_area = [cell_area, repmat({area}, 1, data.N)];
                cell_id = [cell_id, data.cell_id];
                spikes = [spikes;data.spikes];
                for i=1:n_trial
                    rasters{1, i} = [rasters{1, i}; data.rasters{1, i}];
                    firing_rates{1, i} = [firing_rates{1, i}; data.firing_rates{1, i}];
                end
            end
            % save data
            check_path(['../GLM_data/', control, subsession]);
            merged_file = ['../GLM_data/', control, subsession, '/raster_', ...
                    control, subsession, '_', int2str(session_idx), '_0.mat'];
            save(merged_file, 'N', "cell_id", "cuetype", "firing_rates", "n_trial",...
                "rasters", "spikes", "session_name_full", "trial_len", "cell_area");
        end
    end
end
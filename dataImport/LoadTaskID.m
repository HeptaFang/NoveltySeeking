% session_types = {'Muscimol', 'Saline', 'SimRec'};
session_types = {'Muscimol_pre', 'Saline_pre', 'SimRec'};
area_names = {'ACC', 'Thalamus', 'VLPFC'};

unique_sessions_all = ...
    {{'10272023', '11012023', '11102023', '11172023', '12012023',...
    '12082023', '12152023', '12292023', '01052024', '01122024'},...
    {'01302024', '02022024', '02092024', '02162024', '02292024'},...
    {'08112023', '08142023', '08152023', '08162023', '08172023'}};% Muscimol, Saline, SimRec

n_session_types = length(session_types);
n_area_names = length(area_names);

for session_type_idx = 1:n_session_types
    for area_name_idx = 1:n_area_names
        session_type = session_types{session_type_idx};
        area_name = area_names{area_name_idx};
        unique_sessions = unique_sessions_all{session_type_idx};
        n_sessions = length(unique_sessions);

        file_name = ['../taskID/SI_', area_name, '_', session_type, '_popPav_01242025.mat'];
        load(file_name, 'pop');
        




    end
end
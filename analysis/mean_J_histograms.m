function mean_J_histograms(session_type, prepost)
root_path = '\\storage1.ris.wustl.edu\ilyamonosov\Active\Tianhong\';
kernel = 'Delta';
reg = 'L2=2';
epoch = '2500';
area_names = {'ACC', 'Thalamus', 'VLPFC'};
% load data
states = {'Decision', 'Info', 'InfoAnti', 'InfoResp', 'RestClose', 'RestOpen'};
n_states = length(states);
if strcmp(session_type, 'Muscimol')
    n_session = 10;
else
    n_session = 5;
end

% load metadata
if strcmp(prepost, 'Pre')
    session_stage_full = [session_type, 'Pre', states{1}, '_full'];
else
    session_stage_full = [session_type, 'Post', states{1}, '_cortex'];
end

% load model
file_path = [root_path, 'GLM_model/', session_stage_full,...
        '/GLM_', session_stage_full, '_', '1', '_',...
        kernel, '_0_', reg, '_', epoch, '.mat'];
load(file_path,"n_conn_kernel", "kernel_len", "N");

% if strcmp(prepost, 'Post')
%     J_data = cell(2, 2, n_conn_kernel, n_states, n_session);
% else
%     J_data = cell(3, 3, n_conn_kernel, n_states, n_session);
% end

all_J = [];
J_state = [];
J_session = [];
J_kernel = [];
J_i = [];
J_j = [];

for session_idx = 1:n_session
    for state_idx = 1:n_states
        state = states{state_idx};
        if strcmp(prepost, 'Pre')
            session_stage_full = [session_type, 'Pre', state, '_full'];
        else
            session_stage_full = [session_type, 'Post', state, '_cortex'];
        end

        file_path = [root_path, 'GLM_model/', session_stage_full,...
            '/GLM_', session_stage_full, '_', num2str(session_idx), '_',...
            kernel, '_0_', reg, '_', epoch, '.mat'];

        load(file_path, "model_par", "n_PS_kernel", "kernel_len", "N");
        load([root_path, 'GLM_data/', session_stage_full,'/borders_', session_stage_full, '_', ...
                num2str(session_idx),'.mat'], "borders");
        
        borders = [1, borders-0.5]; % area i is from borders(i) to borders(i+1)
        n_area = length(borders) - 1;
        if strcmp(prepost, 'Post')
            assert(n_area == 2, 'Only 2 areas in Post sessions');
        else
            assert(n_area == 3, 'Only 3 areas in Pre sessions');
        end

        for kernel_idx = 1:n_conn_kernel
            J_mat = model_par(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));
            for i = 1:n_area
                for j = 1:n_area
                    J_area = J_mat(borders(i):borders(i+1)-1, borders(j):borders(j+1)-1);
                    % J_data{i, j, kernel_idx, state_idx, session_idx} = J_area(:);
                    all_J = [all_J; J_area(:)];
                    J_state = [J_state; repelem(state_idx, numel(J_area))];
                    J_session = [J_session; repelem(session_idx, numel(J_area))];
                    J_kernel = [J_kernel; repelem(kernel_idx, numel(J_area))];
                    J_i = [J_i; repelem(i, numel(J_area))];
                    J_j = [J_j; repelem(j, numel(J_area))];
                end
            end
        end
    end
end

if strcmp(prepost, 'Pre')
    session_stage_full = [session_type, 'Pre', '_full'];
else
    session_stage_full = [session_type, 'Post', '_cortex'];
end

for kernel_idx = 1:n_conn_kernel
    % Each kernel is one figure
    f = figure("Visible", "off");
    % ij are subplots
    tiledlayout(n_area, n_area);

    for i = 1:n_area
        for j = 1:n_area
            nexttile;
            % violin plot
            filter = J_i == i & J_j == j & J_kernel == kernel_idx;
            data = all_J(filter);
            data_state = J_state(filter);
            data_session = J_session(filter);
            violin(data_session, data, GroupByColor=data_state);
            legend(states);
            title([area_names{j}, ' to ', area_names{i}]);
        end
    end
    suptitle(['Kernel ', num2str(kernel_idx)]);

    fig_folder = [root_path, 'figures/GLM/', session_stage_full];
    check_path(fig_folder);
    saveas(f, [fig_folder, '/J_histograms_kernel_', num2str(kernel_idx), '.png']);
end
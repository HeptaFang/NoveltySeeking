session_type = 'Muscimol';
% session_type = 'Saline';
root_path = '../';
kernel = 'DeltaPure';
reg = 'L2=2';
epoch = '2500';
area_names = {'ACC', 'Thalamus', 'VLPFC'};
% load data
% states = {'Offer1', 'Offer2', 'Decision', 'InfoAnti', 'InfoResp', 'Reward', 'RandomA', 'RandomB'};
% states = {'Offer1', 'Offer2', 'Decision', 'RandomShort', 'RandomLong', 'RandomA', 'RandomB', 'RestOpen', 'RestClose'};
% states = {'RandomA', 'RandomShort', 'RandomLong'};
% states = {'Simulated', 'Simulated_higher', 'RandomLong'};
states = {'Task', 'RestOpen', 'RestClose'};
aligns = {'AlignLast', 'AlignFirst', 'AlignRandom'};
n_states = length(states);
area_mode = 'cortex';

filter_threshold = 1; % J is significant if |J|>filter_threshold*error

if strcmp(session_type, 'Muscimol')
    n_session = 10;
else
    n_session = 5;
end

n_conn_kernel = 3;
J_data = cell(3, 3, n_conn_kernel, n_states, 3, n_session, 2); % (area i, area j, kernel, state, align, session, prepost), each cell is a n_area_i x n_area_j matrix
J_data_err = cell(3, 3, n_conn_kernel, n_states, 3, n_session, 2);

all_J = [];
J_state = [];
J_session = [];
J_kernel = [];
J_i = [];
J_j = [];

prepost_all = {'Pre', 'Post', 'Pre_partial', 'Pre_cortex'};
for session_idx = 1:n_session
    fprintf('Loading session %d\n', session_idx);
    for state_idx = 1:n_states
        state = states{state_idx};
        
        for kernel_idx = 1:n_conn_kernel

            f = figure('Position', [100, 100, 1600, 1600], "Visible", "off");
            t = tiledlayout(2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

            for prepost_idx = 1:4
                prepost = prepost_all{prepost_idx};
                if strcmp(prepost, 'Pre') || strcmp(prepost, 'Pre_partial')
                    session_stage = [session_type, 'Pre', state, '_full'];
                    % session_stage = [session_type, 'Pre', state, '_cortex'];
                elseif strcmp(prepost, 'Post')
                    session_stage = [session_type, 'Post', state, '_cortex'];
                else
                    session_stage = [session_type, 'Pre', state, '_cortex'];
                end
                
                ax = nexttile(t, prepost_idx);
                for align_idx = 1:1
                    align = aligns{align_idx};
                    session_stage_full = [session_stage, '_', align];
                    file_path = [root_path, 'GLM_model/', session_stage_full,...
                        '/GLM_', session_stage_full, '_', num2str(session_idx), '_',...
                        kernel, '_0_', reg, '_', epoch, '.mat'];
            
                    fprintf('Loading %s\n', session_stage_full);
            
                    load(file_path, "model_par", "n_PS_kernel", "kernel_len", "N", "model_err");
                    load([root_path, 'GLM_data/', session_stage_full,'/borders_', session_stage_full, '_', ...
                            num2str(session_idx),'.mat'], "borders");

                    % load sort_idx
                    load([root_path, 'GLM_data/', session_stage_full,'/sortidx_', session_stage_full, '_', ...
                            num2str(session_idx),'_DeltaPure.mat'], "sort_idx");

                    % load cell properties
                    load([root_path, 'GLM_data/', session_stage_full,'/raster_', session_stage_full, '_', ...
                            num2str(session_idx),'_0.mat'], "cell_area", "cell_id", "session_name_full");
                    session_name = session_name_full(1:8);
                    
                    borders = [1, borders+0.5]; % area i is from borders(i) to borders(i+1)
                    if strcmp(prepost, 'Post') || strcmp(prepost, 'Pre_cortex')
                        borders = [borders(1), borders(2), borders(2:end)]; % duplicate the second border for Thalamus
                    end
                    n_area = length(borders) - 1;
                    % if strcmp(prepost, 'Post')
                    %     assert(n_area == 2, 'Only 2 areas in Post sessions');
                    % else
                    %     assert(n_area == 3, 'Only 3 areas in Pre sessions');
                    % end
            
                    J_mat = model_par(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));
                    J_mat = J_mat(sort_idx, sort_idx);
                    % J_err = model_err.total(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));

                    % temp use: fix this
                    if isa(model_err, 'struct')
                        % model_err = model_err.minuslogL;
                        model_err = model_err.total;
                    end

                    J_err = model_err(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));
                    J_err = J_err(sort_idx, sort_idx);
                    omitThr = strcmp(prepost, 'Pre_partial');
                    % Set node_colors
                    node_colors = zeros(N, 3);

                    % load encoding
                    cell_encodings = cell(1, 3);
                    if strcmp(prepost, 'Post')
                        encoding_prepost = 'post';
                        encoding_prepost = 'pre';
                    else
                        encoding_prepost = 'pre';
                    end

                    for area_idx = 1:3
                        area_name = area_names{area_idx};
                        % file_path = [root_path, 'fromMengxi/Sl_', area_name, '_', session_type, '_prevspost_unitSelect.mat'];
                        % load(file_path, 'uSelect');
                        % cell_encodings{area_idx} = uSelect;
                        file_path = [root_path, 'fromMengxi/Sl_', area_name, '_', session_type, '_', encoding_prepost, '_popPav_02202025.mat'];
                        load(file_path, 'pop');
                        cell_encodings{area_idx} = pop;
                    end

                    for i = 1:N
                        if i < borders(2)
                            encoding = cell_encodings{1};
                        elseif i < borders(3)
                            encoding = cell_encodings{2};
                        else
                            encoding = cell_encodings{3};
                        end
                        cell_name = [session_name, '-008_', cell_id{i}];

                        % find idx in encoding.unitID
                        idx = find(strcmp(encoding.unitID, cell_name));
                        if isempty(idx)
                            node_colors(i, :) = [0.25, 0.25, 0.25];
                        else
                            % use the encoding you want
                            if strcmp(encoding.ID.IvsCer_sign{idx}, 'Posi')
                                node_colors(i, :) = [1, 0, 0];
                            elseif strcmp(encoding.ID.IvsCer_sign{idx}, 'Nega')
                                node_colors(i, :) = [0, 0, 1];
                            else
                                node_colors(i, :) = [1, 1, 1];
                            end
                            end
                    end

                    network_plot_hemi(ax, J_mat, J_err, borders, node_colors, omitThr, area_mode, 'full');
                end
                % title(ax, [session_stage_full, ' ', prepost]);
                title(ax, [prepost]);
            end
            title(t, [state, ' session ', num2str(session_idx), ' kernel ', num2str(kernel_idx)]);
            saveas(f, [root_path, 'figures/GLM/network/', session_stage, '_', area_mode, '_', num2str(session_idx), '_', num2str(kernel_idx), '.png']);
        end
    end
end
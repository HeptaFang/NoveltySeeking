% Load all J data and save to a cell array
% J_data: (area i, area j, kernel, state, align, session, prepost), each cell is a n_area_i x n_area_j matrix
% J_data_err: (area i, area j, kernel, state, align, session, prepost), each cell is a n_area_i x n_area_j matrix
% J_cell: (area i, in/out/mean/concatenated, kernel, state, align, session, prepost), each cell is a n_area_i x n_area_j matrix
% J_cell_err: (area i, in/out/mean/concatenated, kernel, state, align, session, prepost), each cell is a n_area_i x n_area_j matrix
% firing_rate: (area, state, align, session, prepost), each cell is a n_area_i x 1 vector

%% parameters
% session_type = 'Muscimol';
session_type = 'Saline';
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
% aligns = {'AlignFirst', 'AlignLast', 'AlignRandom'};
aligns = {'AlignLast'};
n_states = length(states);
n_aligns = length(aligns);
if strcmp(session_type, 'Muscimol')
    n_session = 10;
else
    n_session = 5;
end

%% Load J matrix and error
n_conn_kernel = 3;
n_areas = 3;
J_area = cell(n_areas, n_areas, n_conn_kernel, n_states, n_aligns, n_session, 3); % (area i, area j, kernel, state, align, session, prepost), each cell is a n_area_i x n_area_j matrix
J_area_err = cell(n_areas, n_areas, n_conn_kernel, n_states, n_aligns, n_session, 3);
J_cell = cell(n_areas, 4, n_conn_kernel, n_states, n_aligns, n_session, 3);% (area i, in/out/mean/concatenated, kernel, state, align, session, prepost), each cell is a n_area_i x n_area_j matrix
J_cell_err = cell(n_areas, 4, n_conn_kernel, n_states, n_aligns, n_session, 3);
firing_rate = cell(n_areas, n_states, n_aligns, n_session, 3); % (area, state, align, session, prepost), each cell is a n_area_i x 1 vector

all_J = [];
J_state = [];
J_session = [];
J_kernel = [];
J_i = [];
J_j = [];

prepost_all = {'Pre', 'Post', 'Pre_full'}; % pre cortex, post cortex, pre full
for prepost_idx = 1:3
    prepost = prepost_all{prepost_idx};
    for session_idx = 1:n_session
        fprintf('Loading session %d\n', session_idx);
        for state_idx = 1:n_states
            state = states{state_idx};
            if strcmp(prepost, 'Pre')
                % session_stage = [session_type, 'Pre', state, '_full'];
                session_stage = [session_type, 'Pre', state, '_cortex'];
            elseif strcmp(prepost, 'Post')
                session_stage = [session_type, 'Post', state, '_cortex'];
            else
                session_stage = [session_type, 'Pre', state, '_full'];
            end

            for align_idx = 1:n_aligns
                align = aligns{align_idx};
                session_stage_full = [session_stage, '_', align];
                file_path = [root_path, 'GLM_model/', session_stage_full,...
                    '/GLM_', session_stage_full, '_', num2str(session_idx), '_',...
                    kernel, '_0_', reg, '_', epoch, '.mat'];
        
                fprintf('Loading %s\n', session_stage_full);
        
                load(file_path, "model_par", "n_PS_kernel", "kernel_len", "N", "model_err");
                load([root_path, 'GLM_data/', session_stage_full,'/borders_', session_stage_full, '_', ...
                        num2str(session_idx),'.mat'], "borders");
                
                borders = [1, borders+0.5]; % area i is from borders(i) to borders(i+1)
                n_area = length(borders) - 1;
                % if strcmp(prepost, 'Post')
                %     assert(n_area == 2, 'Only 2 areas in Post sessions');
                % else
                %     assert(n_area == 3, 'Only 3 areas in Pre sessions');
                % end

                % firing rate
                file_path = [root_path, 'GLM_data/', session_stage_full,'/GLMdata_', session_stage_full, ...
                    '_', num2str(session_idx), '_', kernel, '_0.mat'];
                load(file_path, 'raster');
                for i = 1:n_area
                    if n_area == 2 && i == 2
                        % if i == 2
                            i_eff = 3;
                        else
                            i_eff = i;
                    end
                    firing_rate{i_eff, state_idx, align_idx, session_idx, prepost_idx} = mean(raster(borders(i):borders(i+1)-1, :), 2)*1000;
                end
        
                for kernel_idx = 1:n_conn_kernel
                    J_all = model_par(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));
                    % J_err = model_err.total(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));

                    % temp use: fix this
                    if isa(model_err, 'struct')
                        % model_err = model_err.minuslogL;
                        model_err = model_err.total;
                    end

                    J_err_all = model_err(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));

                    for i = 1:n_area
                        for j = 1:n_area
                            if n_area == 2 && i == 2
                            % if i == 2
                                i_eff = 3;
                            else
                                i_eff = i;
                            end
                            if n_area == 2 && j == 2
                            % if j == 2
                                j_eff = 3;
                            else
                                j_eff = j;
                            end

                            J = J_all(borders(i):borders(i+1)-1, borders(j):borders(j+1)-1);
                            J_err = J_err_all(borders(i):borders(i+1)-1, borders(j):borders(j+1)-1);
                            J_area{i_eff, j_eff, kernel_idx, state_idx, align_idx, session_idx, prepost_idx} = J;
                            J_area_err{i_eff, j_eff, kernel_idx, state_idx, align_idx, session_idx, prepost_idx} = J_err;
                            % all_J = [all_J; J_area(:)];
                            % J_state = [J_state, repelem(state_idx, numel(J_area))];
                            % J_session = [J_session, repelem(session_idx, numel(J_area))];
                            % J_kernel = [J_kernel, repelem(kernel_idx, numel(J_area))];
                            % J_i = [J_i, repelem(i_eff, numel(J_area))];
                            % J_j = [J_j, repelem(j_eff, numel(J_area))];
                        end
                    end

                    for i = 1:n_area
                        if n_area == 2 && i == 2
                            % if i == 2
                                i_eff = 3;
                            else
                                i_eff = i;
                        end
                        for direction = 1:2 % in or out
                            if direction == 1 % in
                                J = J_all(borders(i):borders(i+1)-1, :);
                                J_err = J_err_all(borders(i):borders(i+1)-1, :);
                            else % out
                                J = J_all(:, borders(i):borders(i+1)-1).';
                                J_err = J_err_all(:, borders(i):borders(i+1)-1).';
                            end
                            J_cell{i_eff, direction, kernel_idx, state_idx, align_idx, session_idx, prepost_idx} = J;
                            J_cell_err{i_eff, direction, kernel_idx, state_idx, align_idx, session_idx, prepost_idx} = J_err;
                        end
                        J_cell{i_eff, 3, kernel_idx, state_idx, align_idx, session_idx, prepost_idx} = (J_cell{i_eff, 1, kernel_idx, state_idx, align_idx, session_idx, prepost_idx} + J_cell{i_eff, 2, kernel_idx, state_idx, align_idx, session_idx, prepost_idx})/2;
                        J_cell_err{i_eff, 3, kernel_idx, state_idx, align_idx, session_idx, prepost_idx} = sqrt(J_cell_err{i_eff, 1, kernel_idx, state_idx, align_idx, session_idx, prepost_idx}.^2 + J_cell_err{i_eff, 2, kernel_idx, state_idx, align_idx, session_idx, prepost_idx}.^2)/2;
                        J_cell{i_eff, 4, kernel_idx, state_idx, align_idx, session_idx, prepost_idx} = [J_cell{i_eff, 1, kernel_idx, state_idx, align_idx, session_idx, prepost_idx}, J_cell{i_eff, 2, kernel_idx, state_idx, align_idx, session_idx, prepost_idx}];
                        J_cell_err{i_eff, 4, kernel_idx, state_idx, align_idx, session_idx, prepost_idx} = [J_cell_err{i_eff, 1, kernel_idx, state_idx, align_idx, session_idx, prepost_idx}, J_cell_err{i_eff, 2, kernel_idx, state_idx, align_idx, session_idx, prepost_idx}];
                    end
                end
            end
        end
    end
end

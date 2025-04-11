% Use last aligned time bins
% Complete figure of connectivity change
% pre-post connection count of Saline and Muscimol, and relative change
% three groups: ACC, VLPFC, Across
% output = (3 states x 3 kernels) figures. 2x2 subplots.

root_path = '../';
kernel = 'DeltaPure';
reg = 'L2=2';
epoch = '2500';
area_names = {'ACC', 'Thalamus', 'VLPFC'};
area_type_names = {'ACC', 'VLPFC', 'Across area', 'Within area'};
filter_threshold = 1;

% load data
states = {'Task', 'RestOpen', 'RestClose'};
aligns = {'AlignLast'};
session_types = {'Muscimol', 'Saline'};
n_states = length(states);
n_session = 10;
n_conn_kernel = 3;

% (area i, area j, kernel, state, align, session, prepost, session_type), each cell is a n_area_i x n_area_j matrix
J_data = cell(3, 3, n_conn_kernel, n_states, 1, n_session, 2, 2); J_data_err = cell(3, 3, n_conn_kernel, n_states, 1, n_session, 2, 2);
area_size = zeros(3, 3, n_session);

all_J = [];
J_state = [];
J_session = [];
J_kernel = [];
J_i = [];
J_j = [];
for session_type_idx = 1:2
    session_type = session_types{session_type_idx};
    prepost_all = {'Pre', 'Post'};
    if strcmp(session_type, 'Muscimol')
        n_session = 10;
    else
        n_session = 5;
    end
    for prepost_idx = 1:2
        prepost = prepost_all{prepost_idx};
        for session_idx = 1:n_session
            fprintf('Loading session %d\n', session_idx);
            for state_idx = 1:n_states
                state = states{state_idx};
                if strcmp(prepost, 'Pre')
                    session_stage = [session_type, 'Pre', state, '_full'];
                    % session_stage_full = [session_type, 'Pre', state, '_cortex'];
                else
                    session_stage = [session_type, 'Post', state, '_cortex'];
                end

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
                    
                    borders = [1, borders+0.5]; % area i is from borders(i) to borders(i+1)
                    n_area = length(borders) - 1;
            
                    for kernel_idx = 1:n_conn_kernel
                        J_mat = model_par(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));

                        % temp use: fix this
                        if isa(model_err, 'struct')
                            % model_err = model_err.minuslogL;
                            model_err = model_err.total;
                        end

                        J_err = model_err(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));
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

                                J_area = J_mat(borders(i):borders(i+1)-1, borders(j):borders(j+1)-1);
                                J_area_err = J_err(borders(i):borders(i+1)-1, borders(j):borders(j+1)-1);
                                J_data{i_eff, j_eff, kernel_idx, state_idx, align_idx, session_idx, prepost_idx, session_type_idx} = J_area;
                                J_data_err{i_eff, j_eff, kernel_idx, state_idx, align_idx, session_idx, prepost_idx, session_type_idx} = J_area_err;
                                area_size(i_eff, j_eff) = numel(J_area);
                                all_J = [all_J; J_area(:)];
                                J_state = [J_state, repelem(state_idx, numel(J_area))];
                                J_session = [J_session, repelem(session_idx, numel(J_area))];
                                J_kernel = [J_kernel, repelem(kernel_idx, numel(J_area))];
                                J_i = [J_i, repelem(i_eff, numel(J_area))];
                                J_j = [J_j, repelem(j_eff, numel(J_area))];
                            end
                        end
                    end
                end
            end
        end
    end
end

n_area = 3;
ylim_all_kernel = {[-0.05, 0.9], [-0.15, 1], [-0.01, 0.12]};
ylim_all_type = {[], [], [0, 0.075], [0, 0.2]};

for session_type_idx = 1:2
    session_type = session_types{session_type_idx};
    for kernel_idx = 1:n_conn_kernel
        f = figure("Position", [100, 100, 800, 1200], "Visible", "off");
        t = tiledlayout(3, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
        for state_idx = 1:n_states
            state = states{state_idx};
            for area_type = 3:4 % ACC, VLPFC, across area, within_area
                if strcmp(session_type, 'Muscimol')
                    n_session = 10;
                else
                    n_session = 5;
                end

                %% histogram of significant count
                nexttile(t, area_type-2 + (state_idx - 1) * 2);

                % barplot
                sig_ratio = zeros(2, 2); % pos/neg, pre/post
                sig_error = zeros(2, 2); 

                for prepost_idx = 1:2
                    prepost = prepost_all{prepost_idx};

                    total_connection = 0;
                    significant_connection_pos = 0;
                    significant_connection_neg = 0;
                    for session_idx = 1:n_session
                        if area_type == 1 % ACC
                            J_session = J_data{1, 1, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            J_session_err = J_data_err{1, 1, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            conn_session = numel(J_session);
                            sig_session_pos = sum(J_session(:) > J_session_err(:)*filter_threshold);
                            sig_session_neg = sum(J_session(:) < -J_session_err(:)*filter_threshold);
                        elseif area_type == 2 % VLPFC
                            J_session = J_data{3, 3, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            J_session_err = J_data_err{3, 3, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            conn_session = numel(J_session);
                            sig_session_pos = sum(J_session(:) > J_session_err(:)*filter_threshold);
                            sig_session_neg = sum(J_session(:) < -J_session_err(:)*filter_threshold);
                        elseif area_type == 3 % across area
                            J_session = J_data{1, 3, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            J_session_err = J_data_err{1, 3, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            conn_session = numel(J_session);
                            sig_session_pos = sum(J_session(:) > J_session_err(:)*filter_threshold);
                            sig_session_neg = sum(J_session(:) < -J_session_err(:)*filter_threshold);
                            J_session = J_data{3, 1, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            J_session_err = J_data_err{3, 1, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            conn_session = conn_session + numel(J_session);
                            sig_session_pos = sig_session_pos + sum(J_session(:) > J_session_err(:)*filter_threshold);
                            sig_session_neg = sig_session_neg + sum(J_session(:) < -J_session_err(:)*filter_threshold);
                        elseif area_type == 4 % within area
                            J_session = J_data{1, 1, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            J_session_err = J_data_err{1, 1, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            conn_session = numel(J_session);
                            sig_session_pos = sum(J_session(:) > J_session_err(:)*filter_threshold);
                            sig_session_neg = sum(J_session(:) < -J_session_err(:)*filter_threshold);
                            J_session = J_data{3, 3, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            J_session_err = J_data_err{3, 3, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            conn_session = conn_session + numel(J_session);
                            sig_session_pos = sig_session_pos + sum(J_session(:) > J_session_err(:)*filter_threshold);
                            sig_session_neg = sig_session_neg + sum(J_session(:) < -J_session_err(:)*filter_threshold);
                        end
                        total_connection = total_connection + conn_session;
                        significant_connection_pos = significant_connection_pos + sig_session_pos;
                        significant_connection_neg = significant_connection_neg + sig_session_neg;
                    end
                    sig_ratio(1, prepost_idx) = significant_connection_pos / total_connection; % positive connection
                    sig_ratio(2, prepost_idx) = significant_connection_neg / total_connection; % negative connection
                    sig_error(1, prepost_idx) = sqrt(sig_ratio(1, prepost_idx)*(1-sig_ratio(1, prepost_idx))/total_connection);
                    sig_error(2, prepost_idx) = sqrt(sig_ratio(2, prepost_idx)*(1-sig_ratio(2, prepost_idx))/total_connection);
                end

                % barplot
                bar_width = 0.2;
                colors = [1,0,0; 0,0,1]; % pos: red, neg: blue
                lighter_colors = 1-(1-colors)*0.25;

                hold on;
                for prepost_idx = 1:2
                    for posneg = 1:2
                        x = prepost_idx + (posneg-1.5)*bar_width*1.2;
                        bar(x, sig_ratio(posneg, prepost_idx), 'FaceColor', colors(posneg, :), 'BarWidth', bar_width); 
                        errorbar(x, sig_ratio(posneg, prepost_idx), sig_error(posneg, prepost_idx), 'k', 'LineStyle', 'none', 'LineWidth', 1);
                    end
                end
                hold off;

                % set axis limits and labels
                xlim([0.5, 2.5]);
                xticks(1:2);
                xticklabels({'Pre', 'Post'});
                title([state, ', ', area_type_names{area_type}]);
                % legend({'Positive', '', 'Negative'}, 'Location', 'northwest');
                ylabel('Significant connection Ratio');
                set(gca, 'FontSize', 12);
                % set(gca, 'XTickLabelRotation', 45);
                max_y = max(sig_ratio(:) + sig_error(:)) * 1.2;
                ylim([0, max_y]);
                ylim(ylim_all_type{area_type});

            end
        end
        title(t, [session_type, ', Kernel ', num2str(kernel_idx)]);
        fig_folder = [root_path, 'figures/GLM/J_Count'];
        check_path(fig_folder);
        saveas(f, [fig_folder, '/J_Count_', session_type, '_', num2str(kernel_idx), '.png']);
    end
end
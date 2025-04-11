% Use last aligned time bins
% Complete figure of connectivity change
% pre-post connection count of Saline and Muscimol, and relative change
% three groups: ACC, VLPFC, Across
% use all data to one figure: relative sig count post/pre

root_path = '../';
kernel = 'DeltaPure';
reg = 'L2=2';
epoch = '2500';
area_names = {'ACC', 'Thalamus', 'VLPFC'};
area_type_names = {'Within ACC', 'Within VLPFC', 'Across area', 'Within ACC/VLPFC'};
filter_threshold = 1;

% load data
states = {'All', 'Task', 'RestOpen', 'RestClose'};
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
ylim_all_type = {[], [], [0, 0.04], [0, 0.15]};

for kernel_idx = 1:n_conn_kernel
    
    f = figure("Position", [100, 100, 900, 1200], "Visible", "off");
    t = tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    for area_type_idx = 3:4 % ACC, VLPFC, across area, within_area
        nexttile(t, area_type_idx-2);
        plot_data = zeros(2, n_states, 2); % muscimol/saline, state, pos/neg
        plot_data_err = zeros(2, n_states, 2); % muscimol/saline, state, pos/neg

        p_data = zeros(2, n_states, 2); % muscimol/saline, state, pos/neg
        CI_data = cell(2, n_states, 2); % bootstrap CI, muscimol/saline, state, pos/neg
        bootstrap_p = zeros(n_states, 2); % bootstrap p, comparasion of muscimol vs saline, (state, pos/neg)

        for state_idx = 1:n_states
            state = states{state_idx};
            N_all = zeros(2, 2, 2); % total connection num, muscimol/saline, pre/post, pos/neg
            M_all = zeros(2, 2, 2); % significant connection num, muscimol/saline, pre/post, pos/neg

            for session_type_idx = 1:2
                session_type = session_types{session_type_idx};
                if strcmp(session_type, 'Muscimol')
                    n_session = 10;
                else
                    n_session = 5;
                end
                % barplot
                sig_ratio = zeros(2, 2); % pos/neg, pre/post
                sig_error = zeros(2, 2); 
                % sig_diff = zeros(2, 2); % pos/neg
                % sig_diff_error = zeros(2, 2); % pos/neg
                significant_connection = cell(2, 2); % pre/post, pos/neg

                for prepost_idx = 1:2
                    prepost = prepost_all{prepost_idx};

                    for session_idx = 1:n_session
                        if area_type_idx == 1 % ACC
                            J_session = J_data{1, 1, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            J_session_err = J_data_err{1, 1, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            conn_session = numel(J_session);
                            sig_session_pos = J_session(:) > J_session_err(:)*filter_threshold;
                            sig_session_neg = J_session(:) < -J_session_err(:)*filter_threshold;
                        elseif area_type_idx == 2 % VLPFC
                            J_session = J_data{3, 3, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            J_session_err = J_data_err{3, 3, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            conn_session = numel(J_session);
                            sig_session_pos = J_session(:) > J_session_err(:)*filter_threshold;
                            sig_session_neg = J_session(:) < -J_session_err(:)*filter_threshold;
                        elseif area_type_idx == 3 % across area
                            J_session = J_data{1, 3, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            J_session_err = J_data_err{1, 3, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            conn_session = numel(J_session);
                            sig_session_pos = J_session(:) > J_session_err(:)*filter_threshold;
                            sig_session_neg = J_session(:) < -J_session_err(:)*filter_threshold;
                            J_session = J_data{3, 1, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            J_session_err = J_data_err{3, 1, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            conn_session = conn_session + numel(J_session);
                            sig_session_pos = [sig_session_pos; J_session(:) > J_session_err(:)*filter_threshold];
                            sig_session_neg = [sig_session_neg; J_session(:) < -J_session_err(:)*filter_threshold];
                        elseif area_type_idx == 4 % within area
                            J_session = J_data{1, 1, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            J_session_err = J_data_err{1, 1, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            conn_session = numel(J_session);
                            sig_session_pos = J_session(:) > J_session_err(:)*filter_threshold;
                            sig_session_neg = J_session(:) < -J_session_err(:)*filter_threshold;
                            J_session = J_data{3, 3, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            J_session_err = J_data_err{3, 3, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx};
                            conn_session = conn_session + numel(J_session);
                            sig_session_pos = [sig_session_pos; J_session(:) > J_session_err(:)*filter_threshold];
                            sig_session_neg = [sig_session_neg; J_session(:) < -J_session_err(:)*filter_threshold];
                        end

                        significant_connection{prepost_idx, 1} = [significant_connection{prepost_idx, 1}; sig_session_pos];
                        significant_connection{prepost_idx, 2} = [significant_connection{prepost_idx, 2}; sig_session_neg];
                    end
                    total_connection_count = numel(significant_connection{prepost_idx, 1});
                    significant_connection_pos = sum(significant_connection{prepost_idx, 1});
                    significant_connection_neg = sum(significant_connection{prepost_idx, 2});
                    sig_ratio(1, prepost_idx) = significant_connection_pos / total_connection_count; % positive connection
                    sig_ratio(2, prepost_idx) = significant_connection_neg / total_connection_count; % negative connection
                    sig_error(1, prepost_idx) = sqrt(sig_ratio(1, prepost_idx)*(1-sig_ratio(1, prepost_idx))/total_connection_count);
                    sig_error(2, prepost_idx) = sqrt(sig_ratio(2, prepost_idx)*(1-sig_ratio(2, prepost_idx))/total_connection_count);
                    N_all(session_type_idx, prepost_idx, 1) = total_connection_count;
                    M_all(session_type_idx, prepost_idx, 1) = significant_connection_pos;
                    N_all(session_type_idx, prepost_idx, 2) = total_connection_count;
                    M_all(session_type_idx, prepost_idx, 2) = significant_connection_neg;
                end

                % % --mcnemar's test
                % a01_pos = sum(significant_connection{1, 1} & ~significant_connection{2, 1}); % pre=1, post=0
                % a10_pos = sum(~significant_connection{1, 1} & significant_connection{2, 1}); % pre=0, post=1
                % a01_neg = sum(significant_connection{1, 2} & ~significant_connection{2, 2}); % pre=1, post=0
                % a10_neg = sum(~significant_connection{1, 2} & significant_connection{2, 2}); % pre=0, post=1

                % % McNemar's test for pos/neg
                % chi_pos = (a01_pos - a10_pos)^2 / (a01_pos + a10_pos);
                % chi_neg = (a01_neg - a10_neg)^2 / (a01_neg + a10_neg);
                % p_pos = 1 - chi2cdf(chi_pos, 1);
                % p_neg = 1 - chi2cdf(chi_neg, 1);

                % % --paired t-test for pre/post
                % [~, p_pos] = ttest2(significant_connection{1, 1}, significant_connection{2, 1});
                % [~, p_neg] = ttest2(significant_connection{1, 2}, significant_connection{2, 2});

                % p_data(session_type_idx, state_idx, 1) = p_pos;
                % p_data(session_type_idx, state_idx, 2) = p_neg;

                sig_diff = sig_ratio(:, 2) ./ sig_ratio(:, 1);
                sig_diff_error = sig_diff .* sqrt(sig_error(:, 1).^2 ./ sig_ratio(:, 1).^2 + sig_error(:, 2).^2 ./ sig_ratio(:, 2).^2);
                plot_data(session_type_idx, state_idx, :) = sig_diff;
                plot_data_err(session_type_idx, state_idx, :) = sig_diff_error;


            end
            % bootstrap comparasion of muscimol vs saline, (state, pos/neg)

            bootstrap_N = 1000;

            muscimol_ratios = zeros(bootstrap_N, 2); % pos/neg
            saline_ratios = zeros(bootstrap_N, 2); % pos/neg

            for i=1:bootstrap_N
                mus_pre_pos = rand((N_all(1, 1, 1)), 1) < (M_all(1, 1, 1)/N_all(1, 1, 1));
                mus_pre_neg = rand((N_all(1, 1, 2)), 1) < (M_all(1, 1, 2)/N_all(1, 1, 2));
                mus_post_pos = rand((N_all(1, 2, 1)), 1) < (M_all(1, 2, 1)/N_all(1, 2, 1));
                mus_post_neg = rand((N_all(1, 2, 2)), 1) < (M_all(1, 2, 2)/N_all(1, 2, 2));
                sal_pre_pos = rand((N_all(2, 1, 1)), 1) < (M_all(2, 1, 1)/N_all(2, 1, 1));
                sal_pre_neg = rand((N_all(2, 1, 2)), 1) < (M_all(2, 1, 2)/N_all(2, 1, 2));
                sal_post_pos = rand((N_all(2, 2, 1)), 1) < (M_all(2, 2, 1)/N_all(2, 2, 1));
                sal_post_neg = rand((N_all(2, 2, 2)), 1) < (M_all(2, 2, 2)/N_all(2, 2, 2));
                muscimol_ratios(i, 1) = nnz(mus_post_pos)/nnz(mus_pre_pos);
                muscimol_ratios(i, 2) = nnz(mus_post_neg)/nnz(mus_pre_neg);
                saline_ratios(i, 1) = nnz(sal_post_pos)/nnz(sal_pre_pos);
                saline_ratios(i, 2) = nnz(sal_post_neg)/nnz(sal_pre_neg);
            end

            ratio_diff = muscimol_ratios - saline_ratios;

            % CI
            CI_data{1, state_idx, 1} = prctile(muscimol_ratios(:, 1), [2.5, 97.5]);
            CI_data{1, state_idx, 2} = prctile(muscimol_ratios(:, 2), [2.5, 97.5]);
            CI_data{2, state_idx, 1} = prctile(saline_ratios(:, 1), [2.5, 97.5]);
            CI_data{2, state_idx, 2} = prctile(saline_ratios(:, 2), [2.5, 97.5]);

            % muscimol p
            p_left_pos = nnz(muscimol_ratios(:, 1) < 1) / bootstrap_N; % nnz: number of non-zero elements
            p_left_neg = nnz(muscimol_ratios(:, 2) < 1) / bootstrap_N;
            p_right_pos = nnz(muscimol_ratios(:, 1) > 1) / bootstrap_N;
            p_right_neg = nnz(muscimol_ratios(:, 2) > 1) / bootstrap_N;
            p_data(1, state_idx, 1) = min(p_left_pos, p_right_pos)*2; % two-tailed test
            p_data(1, state_idx, 2) = min(p_left_neg, p_right_neg)*2;

            % saline p
            p_left_pos = nnz(saline_ratios(:, 1) < 1) / bootstrap_N; % nnz: number of non-zero elements
            p_left_neg = nnz(saline_ratios(:, 2) < 1) / bootstrap_N;
            p_right_pos = nnz(saline_ratios(:, 1) > 1) / bootstrap_N;
            p_right_neg = nnz(saline_ratios(:, 2) > 1) / bootstrap_N;
            p_data(2, state_idx, 1) = min(p_left_pos, p_right_pos)*2;
            p_data(2, state_idx, 2) = min(p_left_neg, p_right_neg)*2;

            % bootstrap p
            p_left_pos = nnz(ratio_diff(:, 1) < 0) / bootstrap_N; % nnz: number of non-zero elements
            p_left_neg = nnz(ratio_diff(:, 2) < 0) / bootstrap_N;
            p_right_pos = nnz(ratio_diff(:, 1) > 0) / bootstrap_N;
            p_right_neg = nnz(ratio_diff(:, 2) > 0) / bootstrap_N;
            bootstrap_p(state_idx, 1) = min(p_left_pos, p_right_pos)*2;
            bootstrap_p(state_idx, 2) = min(p_left_neg, p_right_neg)*2;
        end
        % plot
        hold on;
        for state_idx = 1:n_states
            for posneg_idx = 1:2
                max_CI = 0;
                for session_type_idx = 1:2
                    x = state_idx + (session_type_idx-1.5)*0.1 + (posneg_idx-1.5)*0.25;
                    y = plot_data(session_type_idx, state_idx, posneg_idx);
                    y_err = plot_data_err(session_type_idx, state_idx, posneg_idx);
                    if posneg_idx == 1
                        color = 'r';
                    else
                        color = 'b';
                    end
                    if session_type_idx == 1
                        alpha = 1;
                    else
                        alpha = 0.3;
                    end

                    b = bar(x, y, 0.1, 'FaceColor', color, 'EdgeColor', 'k', 'FaceAlpha', alpha);

                    % errorbar(x, y, y_err, 'k', 'LineStyle', 'none', 'LineWidth', 1); % propagated standard error

                    % CI
                    CI = CI_data{session_type_idx, state_idx, posneg_idx};
                    CI_center = mean(CI);
                    CI_width = (CI(2) - CI(1)) / 2;
                    errorbar(x, CI_center, CI_width, 'k', 'LineStyle', 'none', 'LineWidth', 1); % CI

                    if CI(2) > max_CI
                        max_CI = CI(2);
                    end

                    % significance for this bar
                    p = p_data(session_type_idx, state_idx, posneg_idx);
                    if p < 0.05
                        if p < 0.001
                            sig_str = '***';
                        elseif p < 0.01
                            sig_str = '**';
                        else
                            sig_str = '*';
                        end
                        % text(x, y + y_err + 0.02, sig_str, 'FontSize', 10, 'Color', 'r', 'HorizontalAlignment', 'center');
                        text(x, CI(2) + 0.02, sig_str, 'FontSize', 10, 'Color', 'r', 'HorizontalAlignment', 'center');
                    
                    end
                end

                % significance check for Muscimol vs Saline
                p = bootstrap_p(state_idx, posneg_idx);
                if p < 0.05
                    if p < 0.001
                        sig_str = '***';
                    elseif p < 0.01
                        sig_str = '**';
                    else
                        sig_str = '*';
                    end
                    x_pos = state_idx + (posneg_idx-1.5)*0.25;
                    y_pos = max_CI + 1;
                    text(x_pos, y_pos, sig_str, 'FontSize', 10, 'Color', 'k', 'HorizontalAlignment', 'center');
                    % plot a line between the two bars
                    x_pos = [x_pos-0.05, x_pos+0.05];
                    y_pos = [y_pos, y_pos];
                    line(x_pos, y_pos, 'Color', 'k', 'LineWidth', 1);
                end

            end

        end
        % dashed line at y=1
        yline(1, 'k--', 'LineWidth', 1);

        hold off;
        xlim([0.5, n_states+0.5]);
        xticks(1:n_states);
        xticklabels(states);
        % xtickangle(45);

        max_y = max(plot_data(:) + plot_data_err(:));
        ylim([0, max_y*1.5]);
        ylabel('significant connection number change (post/pre)');

        title(area_type_names{area_type_idx});
        legends = {'Muscimol positive', '', 'Muscimol negative', '', 'Saline positive', '', 'Saline negative', ''};
        % 2x2 legend, small
        legend(legends, 'Location', 'best', 'NumColumns', 2, 'FontSize', 8, 'Box', 'off', 'Orientation', 'horizontal');
    end
    title(t, ['Kernel ', num2str(kernel_idx)]);
    fig_folder = [root_path, 'figures/GLM/J_Count'];
    check_path(fig_folder);
    saveas(f, [fig_folder, '/J_Count_Relative_', num2str(kernel_idx), '.png']);
end
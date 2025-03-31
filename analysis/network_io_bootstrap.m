% Use last aligned time bins
% plot the bootstrap distribution

session_type = 'Muscimol';
% session_type = 'Saline';
root_path = '../';
kernel = 'DeltaPure';
reg = 'L2=2';
epoch = '2500';
area_names = {'ACC', 'Thalamus', 'VLPFC'};
states = {'Task', 'RestOpen', 'RestClose'};
aligns = {'AlignLast'};
n_states = length(states);
significant_threshold = 1;
use_abs_mean = false;
use_abs_sigmean = false;
% sign_filter = 'raw'; % 'abs', 'pos', 'neg', or 'raw'
plot_legend = true;
plot_stat = true;
bootstrap = true;
flag=3;

if strcmp(session_type, 'Muscimol')
    n_session = 10;
else
    n_session = 5;
end

warning('off');

% sign_filters = {'raw', 'abs', 'pos', 'neg'};
sign_filters = {'raw', 'abs'};
% sign_filters = {'raw'};
for sign_filter_idx = 1:length(sign_filters)
    sign_filter = sign_filters{sign_filter_idx};
    n_conn_kernel = 3;
    J_data = cell(3, 3, n_conn_kernel, n_states, 1, n_session, 2); % (area i, area j, kernel, state, align, session, prepost), each cell is a n_area_i x n_area_j matrix
    J_data_err = cell(3, 3, n_conn_kernel, n_states, 1, n_session, 2);
    area_size = zeros(3, 3, n_session);
    firing_rates = cell(n_states, 1, n_session, 2);

    all_J = [];
    J_state = [];
    J_session = [];
    J_kernel = [];
    J_i = [];
    J_j = [];

    fprintf("%s, Loading data\n", sign_filter);

    prepost_all = {'Pre', 'Post'};
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

                    file_path = [root_path, 'GLM_data/', session_stage_full,...
                        '/GLMdata_', session_stage_full, '_', num2str(session_idx), '_',...
                        kernel, '_0.mat'];
                    load(file_path, "raster")
                    firing_rate = sum(raster, 2) * 1000;
                    firing_rates{state_idx, 1, session_idx, prepost_idx} = firing_rate;
                    
                    borders = [1, borders+0.5]; % area i is from borders(i) to borders(i+1)
                    n_area = length(borders) - 1;
                    % if strcmp(prepost, 'Post')
                    %     assert(n_area == 2, 'Only 2 areas in Post sessions');
                    % else
                    %     assert(n_area == 3, 'Only 3 areas in Pre sessions');
                    % end
            
                    for kernel_idx = 1:n_conn_kernel
                        J_mat = model_par(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));
                        % J_err = model_err.total(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));

                        % temp use: fix this
                        if isa(model_err, 'struct')
                            % model_err = model_err.minuslogL;
                            model_err = model_err.total;
                        end
                        J_err = model_err(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));

                        % set diagonal to NaN
                        J_mat(logical(eye(size(J_mat)))) = NaN;
                        J_err(logical(eye(size(J_err)))) = NaN;

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

                                % sign filter
                                if strcmp(sign_filter, 'abs')
                                    J_area = abs(J_area);
                                elseif strcmp(sign_filter, 'pos')
                                    J_area(J_area < 0) = NaN;
                                elseif strcmp(sign_filter, 'neg')
                                    J_area(J_area > 0) = NaN;
                                end

                                J_data{i_eff, j_eff, kernel_idx, state_idx, align_idx, session_idx, prepost_idx} = J_area;
                                J_data_err{i_eff, j_eff, kernel_idx, state_idx, align_idx, session_idx, prepost_idx} = J_area_err;
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

    % if strcmp(prepost, 'Pre')
    %     session_stage_full = [session_type, 'Pre', '_full'];
    % else
    %     session_stage_full = [session_type, 'Post', '_cortex'];
    % end
    n_area = 3;
    ylim_all_kernel = {[-0.05, 0.9], [-0.15, 1], [-0.01, 0.12]};

    % data_mat = J_data{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx};
    % error_mat = J_data_err{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx};
    for state_idx = 1:n_states
        for align_idx = 1:align_idx
            for kernel_idx = 1:n_conn_kernel
                %% calculate cell profile
                fprintf("%s, %s, Kernel %d, calculating cell profile\n", sign_filter, states{state_idx}, kernel_idx);

                all_profile = cell(1, 0);
                
                for session_idx = 1:n_session
                    N_ACC = size(J_data{1, 1, kernel_idx, 1, 1, session_idx, 1}, 1);
                    N_thalamus = size(J_data{2, 2, kernel_idx, 1, 1, session_idx, 1}, 1);
                    N_VLPFC = size(J_data{3, 3, kernel_idx, 1, 1, session_idx, 1}, 1);
                    N_pre = N_ACC + N_thalamus + N_VLPFC;
                    N_post = N_ACC + N_VLPFC;

                    cell_profile = cell(1, N_pre);
                    area_labels = {'ACC', 'Thalamus', 'VLPFC'};

                    for i = 1:N_pre
                        cell_profile{i}.session = session_idx;
                        if i <= N_ACC
                            cell_profile{i}.area = 1;
                            idx_within_area = i;
                        elseif i <= N_ACC + N_thalamus
                            cell_profile{i}.area = 2;
                            idx_within_area = i - N_ACC;
                        else
                            cell_profile{i}.area = 3;
                            idx_within_area = i - N_ACC - N_thalamus;
                        end

                        % idx: pre/post, target area, in/out/total
                        cell_profile{i}.mean = zeros(2, 3, 3); % mean J to an area
                        cell_profile{i}.sig_count = zeros(2, 3, 3); % number of Sig J to an area
                        cell_profile{i}.sig_mean = zeros(2, 3, 3); % mean of Sig J to an area
                        cell_profile{i}.sig_fixed_mean = zeros(2, 3, 3); % mean J, non-Sig J set to 0
                        cell_profile{i}.all_data = zeros(2, 3, 3, 4); % all J data

                        % calculate mean, significant count, significant mean
                        for prepost_idx = 1:2
                            for target_idx = 1:3
                                if prepost_idx == 2 && (target_idx == 2 || cell_profile{i}.area == 2)
                                    continue;
                                end

                                total_J = 0;
                                total_count = 0;
                                total_sig_J = 0;
                                total_sig_count = 0;

                                % inward connections
                                connections = J_data{cell_profile{i}.area, target_idx, kernel_idx, state_idx, align_idx, session_idx, prepost_idx}(idx_within_area, :);
                                error = J_data_err{cell_profile{i}.area, target_idx, kernel_idx, state_idx, align_idx, session_idx, prepost_idx}(idx_within_area, :);
                                significant_filter = abs(connections) > significant_threshold * error;

                                if use_abs_mean
                                    cell_profile{i}.mean(prepost_idx, target_idx, 1) = mean(abs(connections), "omitmissing");
                                    total_J = total_J + sum(abs(connections), "omitmissing");
                                else
                                    cell_profile{i}.mean(prepost_idx, target_idx, 1) = mean(connections, "omitmissing");
                                    total_J = total_J + sum(connections, "omitmissing");
                                end
                                cell_profile{i}.sig_count(prepost_idx, target_idx, 1) = sum(significant_filter, "omitmissing");
                                if use_abs_sigmean
                                    cell_profile{i}.sig_mean(prepost_idx, target_idx, 1) = mean(abs(connections(significant_filter)), "omitmissing");
                                    total_sig_J = total_sig_J + sum(abs(connections(significant_filter)));
                                else
                                    cell_profile{i}.sig_mean(prepost_idx, target_idx, 1) = mean(connections(significant_filter), "omitmissing");
                                    total_sig_J = total_sig_J + sum(connections(significant_filter), "omitmissing");
                                end

                                total_count = total_count + numel(connections);
                                total_sig_count = total_sig_count + sum(significant_filter);

                                % outward connections
                                connections = J_data{target_idx, cell_profile{i}.area, kernel_idx, state_idx, align_idx, session_idx, prepost_idx}(:, idx_within_area);
                                error = J_data_err{target_idx, cell_profile{i}.area, kernel_idx, state_idx, align_idx, session_idx, prepost_idx}(:, idx_within_area);
                                significant_filter = abs(connections) > significant_threshold * error;

                                if use_abs_mean
                                    cell_profile{i}.mean(prepost_idx, target_idx, 2) = mean(abs(connections), "omitmissing");
                                    total_J = total_J + sum(abs(connections), "omitmissing");
                                else
                                    cell_profile{i}.mean(prepost_idx, target_idx, 2) = mean(connections, "omitmissing");
                                    total_J = total_J + sum(connections, "omitmissing");
                                end
                                cell_profile{i}.sig_count(prepost_idx, target_idx, 2) = sum(significant_filter, "omitmissing");
                                if use_abs_sigmean
                                    cell_profile{i}.sig_mean(prepost_idx, target_idx, 2) = mean(abs(connections(significant_filter)), "omitmissing");
                                    total_sig_J = total_sig_J + sum(abs(connections(significant_filter)), "omitmissing");
                                else
                                    cell_profile{i}.sig_mean(prepost_idx, target_idx, 2) = mean(connections(significant_filter), "omitmissing");
                                    total_sig_J = total_sig_J + sum(connections(significant_filter), "omitmissing");
                                end

                                total_count = total_count + numel(connections);
                                total_sig_count = total_sig_count + sum(significant_filter);

                                % total connections
                                cell_profile{i}.mean(prepost_idx, target_idx, 3) = total_J / total_count;
                                cell_profile{i}.sig_count(prepost_idx, target_idx, 3) = total_sig_count;
                                cell_profile{i}.sig_mean(prepost_idx, target_idx, 3) = total_sig_J / total_sig_count;
                                cell_profile{i}.sig_fixed_mean(prepost_idx, target_idx, 3) = total_sig_J / total_count;

                            end
                        end
                        cell_profile{i}.all_data(:, :, :, 1) = cell_profile{i}.mean(:, :, :);
                        cell_profile{i}.all_data(:, :, :, 2) = cell_profile{i}.sig_count(:, :, :);
                        cell_profile{i}.all_data(:, :, :, 3) = cell_profile{i}.sig_mean(:, :, :);
                        cell_profile{i}.all_data(:, :, :, 4) = cell_profile{i}.sig_fixed_mean(:, :, :);

                        % for cortical neurons, calculate post - pre for within and acrosss
                        if cell_profile{i}.area == 1 || cell_profile{i}.area == 3
                            cell_profile{i}.mean_diff = cell_profile{i}.mean(2, :, :) - cell_profile{i}.mean(1, :, :);
                            cell_profile{i}.sig_count_diff = cell_profile{i}.sig_count(2, :, :) - cell_profile{i}.sig_count(1, :, :);
                            cell_profile{i}.sig_mean_diff = cell_profile{i}.sig_mean(2, :, :) - cell_profile{i}.sig_mean(1, :, :);
                            cell_profile{i}.sig_fixed_mean_diff = cell_profile{i}.sig_fixed_mean(2, :, :) - cell_profile{i}.sig_fixed_mean(1, :, :);
                            cell_profile{i}.all_data_diff = cell_profile{i}.all_data(2, :, :, :) - cell_profile{i}.all_data(1, :, :, :);
                        end
                    end
                    all_profile = [all_profile, cell_profile];
                end

                %% plot
                fprintf("Plotting...");
                % for session_idx = 0:n_session
                for session_idx = 0:0
                    fprintf("%d...", session_idx);
                    %% Figure 1: Pre
                    f = figure("Position", [0, 0, 1800, 1800], "Visible", "off");
                    t = tiledlayout(3, 3, "TileSpacing", "compact", "Padding", "compact");

                    if session_idx == 0
                        cell_profile = all_profile;
                    else
                        filter = cellfun(@(x) x.session == session_idx, all_profile);
                        cell_profile = all_profile(filter);
                    end

                    data_type_labels = {'Mean J', 'Sig J Count', 'Mean of Sig J', 'Mean J, Non-sig set to 0'};

                    %% 1. cortical neurons: pre connection with thalamus vs pre within
                    plot_types = [1, 2, 4]; % 1: mean, 2: significant count, 3: significant mean, 4: significant fixed mean
                    for plot_idx = 1:3
                        nexttile(plot_idx);
                        data_type = plot_types(plot_idx);
                        hold on;
                        % ACC
                        filter = cellfun(@(x) x.area == 1, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_ACC = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                        y_ACC = cellfun(@(x) x.all_data(1, 1, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_ACC) & ~isnan(y_ACC);
                        x_ACC = x_ACC(nan_filter);
                        y_ACC = y_ACC(nan_filter);
                        scatter(x_ACC, y_ACC, 10, 'r', 'x', 'Tag', 'ACC');

                        % VLPFC
                        filter = cellfun(@(x) x.area == 3, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_VLPFC = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                        y_VLPFC = cellfun(@(x) x.all_data(1, 3, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_VLPFC) & ~isnan(y_VLPFC);
                        x_VLPFC = x_VLPFC(nan_filter);
                        y_VLPFC = y_VLPFC(nan_filter);
                        scatter(x_VLPFC, y_VLPFC, 10, 'b', '+', 'Tag', 'VLPFC');

                        x_all = [x_ACC, x_VLPFC];
                        y_all = [y_ACC, y_VLPFC];

                        % linear regression
                        if max(x_all)>min(x_all)
                            x_fit = linspace(min(x_all), max(x_all), 100);
                            plot(x_fit, x_fit, 'k:', 'LineWidth', 1); % identity line
    
                            if length(x_ACC) > 1
                                X = x_ACC;
                                Y = y_ACC;
                                text_pos = [0.5, 0.95];
                                text_color = 'r';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            if length(x_VLPFC) > 1
                                X = x_VLPFC;
                                Y = y_VLPFC;
                                text_pos = [0.5, 0.85];
                                text_color = 'b';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            X = x_all;
                            Y = y_all;
                            text_pos = [0.5, 0.75];
                            text_color = 'm';
                            mdl = fitlm(X, Y);
                            y_fit = predict(mdl, x_fit');
                            plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                            if plot_stat
                                if bootstrap
                                    stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                    r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                    p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                    text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                else
                                    [r, p] = corr(X', Y', 'Type', 'Pearson');
                                    text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                end
                            end
    
                            % text(0.5, 0.85, ['k = ', num2str(mdl_all.Coefficients.Estimate(2))], 'Units', 'normalized');
                            % text(0.5, 0.8, ['b = ', num2str(mdl_all.Coefficients.Estimate(1))], 'Units', 'normalized');
                        end

                        xlabel(['Pre, J to Th, ', data_type_labels{data_type}]);
                        ylabel(['Pre, within, ', data_type_labels{data_type}]);
                        axis equal;

                        % origin axis
                        plot([0, 0], ylim, 'k--');
                        plot(xlim, [0, 0], 'k--');

                        if plot_legend
                            legend({"ACC", "VLPFC"}, "Location", "northwest");
                        end
                    end

                    %% 2. cortical neurons: pre connection with thalamus vs pre across
                    plot_types = [1, 2, 4]; % 1: mean, 2: significant count, 3: significant mean, 4: significant fixed mean
                    for plot_idx = 1:3
                        nexttile(plot_idx+3);
                        data_type = plot_types(plot_idx);
                        hold on;

                        % ACC
                        filter = cellfun(@(x) x.area == 1, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_ACC = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                        y_ACC = cellfun(@(x) x.all_data(1, 3, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_ACC) & ~isnan(y_ACC);
                        x_ACC = x_ACC(nan_filter);
                        y_ACC = y_ACC(nan_filter);
                        scatter(x_ACC, y_ACC, 10, 'r', 'x', 'Tag', 'ACC');

                        % VLPFC
                        filter = cellfun(@(x) x.area == 3, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_VLPFC = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                        y_VLPFC = cellfun(@(x) x.all_data(1, 1, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_VLPFC) & ~isnan(y_VLPFC);
                        x_VLPFC = x_VLPFC(nan_filter);
                        y_VLPFC = y_VLPFC(nan_filter);
                        scatter(x_VLPFC, y_VLPFC, 10, 'b', '+', 'Tag', 'VLPFC');

                        x_all = [x_ACC, x_VLPFC];
                        y_all = [y_ACC, y_VLPFC];

                        % linear regression
                        if max(x_all)>min(x_all)
                            x_fit = linspace(min(x_all), max(x_all), 100);
                            plot(x_fit, x_fit, 'k:', 'LineWidth', 1); % identity line
    
                            if length(x_ACC) > 1
                                X = x_ACC;
                                Y = y_ACC;
                                text_pos = [0.5, 0.95];
                                text_color = 'r';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            if length(x_VLPFC) > 1
                                X = x_VLPFC;
                                Y = y_VLPFC;
                                text_pos = [0.5, 0.85];
                                text_color = 'b';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            X = x_all;
                            Y = y_all;
                            text_pos = [0.5, 0.75];
                            text_color = 'm';
                            mdl = fitlm(X, Y);
                            y_fit = predict(mdl, x_fit');
                            plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                            if plot_stat
                                if bootstrap
                                    stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                    r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                    p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                    text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                else
                                    [r, p] = corr(X', Y', 'Type', 'Pearson');
                                    text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                end
                            end

                            % text(0.5, 0.85, ['k = ', num2str(mdl_all.Coefficients.Estimate(2))], 'Units', 'normalized');
                            % text(0.5, 0.8, ['b = ', num2str(mdl_all.Coefficients.Estimate(1))], 'Units', 'normalized');
                        end

                        xlabel(['Pre, J to Th, ', data_type_labels{data_type}]);
                        ylabel(['Pre, across, ', data_type_labels{data_type}]);
                        axis equal;

                        % origin axis
                        plot([0, 0], ylim, 'k--');
                        plot(xlim, [0, 0], 'k--');

                        if plot_legend
                            legend({"ACC", "VLPFC"}, "Location", "northwest");
                        end
                    end

                    %% 3. cortical neurons: pre within vs pre across
                    plot_types = [1, 2, 4]; % 1: mean, 2: significant count, 3: significant mean, 4: significant fixed mean
                    for plot_idx = 1:3
                        nexttile(plot_idx+6);
                        data_type = plot_types(plot_idx);
                        hold on;

                        % ACC
                        filter = cellfun(@(x) x.area == 1, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_ACC = cellfun(@(x) x.all_data(1, 1, 3, data_type), filtered_profile);
                        y_ACC = cellfun(@(x) x.all_data(1, 3, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_ACC) & ~isnan(y_ACC);
                        x_ACC = x_ACC(nan_filter);
                        y_ACC = y_ACC(nan_filter);
                        scatter(x_ACC, y_ACC, 10, 'r', 'x', 'Tag', 'ACC');

                        % VLPFC
                        filter = cellfun(@(x) x.area == 3, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_VLPFC = cellfun(@(x) x.all_data(1, 3, 3, data_type), filtered_profile);
                        y_VLPFC = cellfun(@(x) x.all_data(1, 1, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_VLPFC) & ~isnan(y_VLPFC);
                        x_VLPFC = x_VLPFC(nan_filter);
                        y_VLPFC = y_VLPFC(nan_filter);
                        scatter(x_VLPFC, y_VLPFC, 10, 'b', '+', 'Tag', 'VLPFC');

                        x_all = [x_ACC, x_VLPFC];
                        y_all = [y_ACC, y_VLPFC];

                        % linear regression
                        if max(x_all)>min(x_all)
                            x_fit = linspace(min(x_all), max(x_all), 100);
                            plot(x_fit, x_fit, 'k:', 'LineWidth', 1); % identity line
    
                            if length(x_ACC) > 1
                                X = x_ACC;
                                Y = y_ACC;
                                text_pos = [0.5, 0.95];
                                text_color = 'r';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            if length(x_VLPFC) > 1
                                X = x_VLPFC;
                                Y = y_VLPFC;
                                text_pos = [0.5, 0.85];
                                text_color = 'b';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            X = x_all;
                            Y = y_all;
                            text_pos = [0.5, 0.75];
                            text_color = 'm';
                            mdl = fitlm(X, Y);
                            y_fit = predict(mdl, x_fit');
                            plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                            if plot_stat
                                if bootstrap
                                    stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                    r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                    p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                    text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                else
                                    [r, p] = corr(X', Y', 'Type', 'Pearson');
                                    text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                end
                            end
                            % text(0.5, 0.85, ['k = ', num2str(mdl_all.Coefficients.Estimate(2))], 'Units', 'normalized');
                            % text(0.5, 0.8, ['b = ', num2str(mdl_all.Coefficients.Estimate(1))], 'Units', 'normalized');
                        end

                        xlabel(['Pre, within, ', data_type_labels{data_type}]);
                        ylabel(['Pre, across, ', data_type_labels{data_type}]);
                        axis equal;

                        % origin axis
                        plot([0, 0], ylim, 'k--');
                        plot(xlim, [0, 0], 'k--');

                        if plot_legend
                            legend({"ACC", "VLPFC"}, "Location", "northwest");
                        end
                    end


                    %% save
                    if session_idx == 0
                        title(t, [session_type, ' All sessions, ', states{state_idx}, ', ', aligns{align_idx}, ', Kernel ', num2str(kernel_idx), ', ', sign_filter]);
                    else
                        title(t, [session_type, 'Session ', num2str(session_idx), ', ', states{state_idx}, ', ', aligns{align_idx}, ', Kernel ', num2str(kernel_idx), ', ', sign_filter]);
                    end

                    save_folder = [root_path, 'figures/network_io/', session_type, '/', sign_filter, '/pre'];
                    if ~exist(save_folder, 'dir')
                        mkdir(save_folder);
                    end
                    save_name = [save_folder, '/network_io_', states{state_idx}, '_', aligns{align_idx}, '_', num2str(session_idx), '_', num2str(kernel_idx), '.png'];
                    if bootstrap && plot_stat
                        save_name = [save_name(1:end-4), '_bootstrap', save_name(end-3:end)];
                    end
                    if ~plot_stat
                        save_name = [save_name(1:end-4), '_no_stat', save_name(end-3:end)];
                    end
                    % saveas(f, save_name);
                    print(f, save_name, '-dpng');
                    
                    close(f);

                    %% Figure 2: Post
                    f = figure("Position", [0, 0, 1800, 1800], "Visible", "off");
                    t = tiledlayout(3, 3, "TileSpacing", "compact", "Padding", "compact");

                    if session_idx == 0
                        cell_profile = all_profile;
                    else
                        filter = cellfun(@(x) x.session == session_idx, all_profile);
                        cell_profile = all_profile(filter);
                    end

                    data_type_labels = {'Mean J', 'Sig J Count', 'Mean of Sig J', 'Mean J, Non-sig set to 0'};

                    %% 1. cortical neurons: pre connection with thalamus vs post within
                    plot_types = [1, 2, 4]; % 1: mean, 2: significant count, 3: significant mean, 4: significant fixed mean
                    for plot_idx = 1:3
                        nexttile(plot_idx);
                        data_type = plot_types(plot_idx);
                        hold on;
                        
                        % ACC
                        filter = cellfun(@(x) x.area == 1, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_ACC = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                        y_ACC = cellfun(@(x) x.all_data(2, 1, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_ACC) & ~isnan(y_ACC);
                        x_ACC = x_ACC(nan_filter);
                        y_ACC = y_ACC(nan_filter);
                        scatter(x_ACC, y_ACC, 10, 'r', 'x', 'Tag', 'ACC');

                        % VLPFC
                        filter = cellfun(@(x) x.area == 3, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_VLPFC = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                        y_VLPFC = cellfun(@(x) x.all_data(2, 3, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_VLPFC) & ~isnan(y_VLPFC);
                        x_VLPFC = x_VLPFC(nan_filter);
                        y_VLPFC = y_VLPFC(nan_filter);
                        scatter(x_VLPFC, y_VLPFC, 10, 'b', '+', 'Tag', 'VLPFC');

                        x_all = [x_ACC, x_VLPFC];
                        y_all = [y_ACC, y_VLPFC];

                        % linear regression
                        if max(x_all)>min(x_all)
                            x_fit = linspace(min(x_all), max(x_all), 100);
                            plot(x_fit, x_fit, 'k:', 'LineWidth', 1); % identity line
    
                            if length(x_ACC) > 1
                                X = x_ACC;
                                Y = y_ACC;
                                text_pos = [0.5, 0.95];
                                text_color = 'r';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            if length(x_VLPFC) > 1
                                X = x_VLPFC;
                                Y = y_VLPFC;
                                text_pos = [0.5, 0.85];
                                text_color = 'b';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            X = x_all;
                            Y = y_all;
                            text_pos = [0.5, 0.75];
                            text_color = 'm';
                            mdl = fitlm(X, Y);
                            y_fit = predict(mdl, x_fit');
                            plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                            if plot_stat
                                if bootstrap
                                    stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                    r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                    p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                    text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                else
                                    [r, p] = corr(X', Y', 'Type', 'Pearson');
                                    text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                end
                            end
    
                            % text(0.5, 0.85, ['k = ', num2str(mdl_all.Coefficients.Estimate(2))], 'Units', 'normalized');
                            % text(0.5, 0.8, ['b = ', num2str(mdl_all.Coefficients.Estimate(1))], 'Units', 'normalized');
                        end

                        xlabel(['Pre, J to Th, ', data_type_labels{data_type}]);
                        ylabel(['Post, within, ', data_type_labels{data_type}]);
                        axis equal;

                        % origin axis
                        plot([0, 0], ylim, 'k--');
                        plot(xlim, [0, 0], 'k--');

                        if plot_legend
                            legend({"ACC", "VLPFC"}, "Location", "northwest");
                        end
                    end

                    %% 2. cortical neurons: pre connection with thalamus vs post across
                    plot_types = [1, 2, 4]; % 1: mean, 2: significant count, 3: significant mean, 4: significant fixed mean
                    for plot_idx = 1:3
                        nexttile(plot_idx+3);
                        data_type = plot_types(plot_idx);
                        hold on;
                        
                        % ACC
                        filter = cellfun(@(x) x.area == 1, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_ACC = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                        y_ACC = cellfun(@(x) x.all_data(2, 3, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_ACC) & ~isnan(y_ACC);
                        x_ACC = x_ACC(nan_filter);
                        y_ACC = y_ACC(nan_filter);
                        scatter(x_ACC, y_ACC, 10, 'r', 'x', 'Tag', 'ACC');

                        % VLPFC
                        filter = cellfun(@(x) x.area == 3, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_VLPFC = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                        y_VLPFC = cellfun(@(x) x.all_data(2, 1, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_VLPFC) & ~isnan(y_VLPFC);
                        x_VLPFC = x_VLPFC(nan_filter);
                        y_VLPFC = y_VLPFC(nan_filter);
                        scatter(x_VLPFC, y_VLPFC, 10, 'b', '+', 'Tag', 'VLPFC');

                        x_all = [x_ACC, x_VLPFC];
                        y_all = [y_ACC, y_VLPFC];

                        % linear regression
                        if max(x_all)>min(x_all)
                            x_fit = linspace(min(x_all), max(x_all), 100);
                            plot(x_fit, x_fit, 'k:', 'LineWidth', 1); % identity line
    
                            if length(x_ACC) > 1
                                X = x_ACC;
                                Y = y_ACC;
                                text_pos = [0.5, 0.95];
                                text_color = 'r';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            if length(x_VLPFC) > 1
                                X = x_VLPFC;
                                Y = y_VLPFC;
                                text_pos = [0.5, 0.85];
                                text_color = 'b';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            X = x_all;
                            Y = y_all;
                            text_pos = [0.5, 0.75];
                            text_color = 'm';
                            mdl = fitlm(X, Y);
                            y_fit = predict(mdl, x_fit');
                            plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                            if plot_stat
                                if bootstrap
                                    stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                    r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                    p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                    text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                else
                                    [r, p] = corr(X', Y', 'Type', 'Pearson');
                                    text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                end
                            end
    
                            % text(0.5, 0.85, ['k = ', num2str(mdl_all.Coefficients.Estimate(2))], 'Units', 'normalized');
                            % text(0.5, 0.8, ['b = ', num2str(mdl_all.Coefficients.Estimate(1))], 'Units', 'normalized');
                        end

                        xlabel(['Pre, J to Th, ', data_type_labels{data_type}]);
                        ylabel(['Post, across, ', data_type_labels{data_type}]);
                        axis equal;

                        % origin axis
                        plot([0, 0], ylim, 'k--');
                        plot(xlim, [0, 0], 'k--');

                        if plot_legend
                            legend({"ACC", "VLPFC"}, "Location", "northwest");
                        end
                    end

                    %% 3. cortical neurons: post within vs post across
                    plot_types = [1, 2, 4]; % 1: mean, 2: significant count, 3: significant mean, 4: significant fixed mean
                    for plot_idx = 1:3
                        nexttile(plot_idx+6);
                        data_type = plot_types(plot_idx);
                        hold on;
                        
                        % ACC
                        filter = cellfun(@(x) x.area == 1, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_ACC = cellfun(@(x) x.all_data(2, 1, 3, data_type), filtered_profile);
                        y_ACC = cellfun(@(x) x.all_data(2, 3, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_ACC) & ~isnan(y_ACC);
                        x_ACC = x_ACC(nan_filter);
                        y_ACC = y_ACC(nan_filter);
                        scatter(x_ACC, y_ACC, 10, 'r', 'x', 'Tag', 'ACC');

                        % VLPFC
                        filter = cellfun(@(x) x.area == 3, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_VLPFC = cellfun(@(x) x.all_data(2, 3, 3, data_type), filtered_profile);
                        y_VLPFC = cellfun(@(x) x.all_data(2, 1, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_VLPFC) & ~isnan(y_VLPFC);
                        x_VLPFC = x_VLPFC(nan_filter);
                        y_VLPFC = y_VLPFC(nan_filter);
                        scatter(x_VLPFC, y_VLPFC, 10, 'b', '+', 'Tag', 'VLPFC');

                        x_all = [x_ACC, x_VLPFC];
                        y_all = [y_ACC, y_VLPFC];

                        % linear regression
                        if max(x_all)>min(x_all)
                            x_fit = linspace(min(x_all), max(x_all), 100);
                            plot(x_fit, x_fit, 'k:', 'LineWidth', 1); % identity line
    
                            if length(x_ACC) > 1
                                X = x_ACC;
                                Y = y_ACC;
                                text_pos = [0.5, 0.95];
                                text_color = 'r';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            if length(x_VLPFC) > 1
                                X = x_VLPFC;
                                Y = y_VLPFC;
                                text_pos = [0.5, 0.85];
                                text_color = 'b';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            X = x_all;
                            Y = y_all;
                            text_pos = [0.5, 0.75];
                            text_color = 'm';
                            mdl = fitlm(X, Y);
                            y_fit = predict(mdl, x_fit');
                            plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                            if plot_stat
                                if bootstrap
                                    stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                    r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                    p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                    text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                else
                                    [r, p] = corr(X', Y', 'Type', 'Pearson');
                                    text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                end
                            end
    
                            % text(0.5, 0.85, ['k = ', num2str(mdl_all.Coefficients.Estimate(2))], 'Units', 'normalized');
                            % text(0.5, 0.8, ['b = ', num2str(mdl_all.Coefficients.Estimate(1))], 'Units', 'normalized');
                        end

                        xlabel(['Post, within, ', data_type_labels{data_type}]);
                        ylabel(['Post, across, ', data_type_labels{data_type}]);
                        axis equal;

                        % origin axis
                        plot([0, 0], ylim, 'k--');
                        plot(xlim, [0, 0], 'k--');

                        if plot_legend
                            legend({"ACC", "VLPFC"}, "Location", "northwest");
                        end
                    end


                    %% save
                    if session_idx == 0
                        title(t, [session_type, ' All sessions, ', states{state_idx}, ', ', aligns{align_idx}, ', Kernel ', num2str(kernel_idx), ', ', sign_filter]);
                    else
                        title(t, [session_type, 'Session ', num2str(session_idx), ', ', states{state_idx}, ', ', aligns{align_idx}, ', Kernel ', num2str(kernel_idx), ', ', sign_filter]);
                    end

                    save_folder = [root_path, 'figures/network_io/', session_type, '/', sign_filter, '/post'];
                    if ~exist(save_folder, 'dir')
                        mkdir(save_folder);
                    end
                    save_name = [save_folder, '/network_io_', states{state_idx}, '_', aligns{align_idx}, '_', num2str(session_idx), '_', num2str(kernel_idx), '.png'];
                    if bootstrap && plot_stat
                        save_name = [save_name(1:end-4), '_bootstrap', save_name(end-3:end)];
                    end
                    if ~plot_stat
                        save_name = [save_name(1:end-4), '_no_stat', save_name(end-3:end)];
                    end
                    % saveas(f, save_name);
                    print(f, save_name, '-dpng');

                    close(f);

                    %% Figure 3: Pre vs Post
                    f = figure("Position", [0, 0, 1200, 1600], "Visible", "off");
                    t = tiledlayout(4, 3, "TileSpacing", "compact", "Padding", "compact");

                    if session_idx == 0
                        cell_profile = all_profile;
                    else
                        filter = cellfun(@(x) x.session == session_idx, all_profile);
                        cell_profile = all_profile(filter);
                    end

                    data_type_labels = {'Mean J', 'Sig J Count', 'Mean of Sig J', 'Mean J, Non-sig set to 0'};

                    %% 1. cortical neurons: pre within vs post within
                    plot_types = [1, 2, 4]; % 1: mean, 2: significant count, 3: significant mean, 4: significant fixed mean
                    for plot_idx = 1:3
                        nexttile(plot_idx);
                        data_type = plot_types(plot_idx);
                        hold on;
                        
                        % ACC
                        filter = cellfun(@(x) x.area == 1, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_ACC = cellfun(@(x) x.all_data(1, 1, 3, data_type), filtered_profile);
                        y_ACC = cellfun(@(x) x.all_data(2, 1, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_ACC) & ~isnan(y_ACC);
                        x_ACC = x_ACC(nan_filter);
                        y_ACC = y_ACC(nan_filter);
                        scatter(x_ACC, y_ACC, 10, 'r', 'x', 'Tag', 'ACC');

                        % VLPFC
                        filter = cellfun(@(x) x.area == 3, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_VLPFC = cellfun(@(x) x.all_data(1, 3, 3, data_type), filtered_profile);
                        y_VLPFC = cellfun(@(x) x.all_data(2, 3, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_VLPFC) & ~isnan(y_VLPFC);
                        x_VLPFC = x_VLPFC(nan_filter);
                        y_VLPFC = y_VLPFC(nan_filter);
                        scatter(x_VLPFC, y_VLPFC, 10, 'b', '+', 'Tag', 'VLPFC');

                        x_all = [x_ACC, x_VLPFC];
                        y_all = [y_ACC, y_VLPFC];

                        % linear regression
                        if max(x_all)>min(x_all)
                            x_fit = linspace(min(x_all), max(x_all), 100);
                            plot(x_fit, x_fit, 'k:', 'LineWidth', 1); % identity line
    
                            if length(x_ACC) > 1
                                X = x_ACC;
                                Y = y_ACC;
                                text_pos = [0.5, 0.95];
                                text_color = 'r';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            if length(x_VLPFC) > 1
                                X = x_VLPFC;
                                Y = y_VLPFC;
                                text_pos = [0.5, 0.85];
                                text_color = 'b';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            X = x_all;
                            Y = y_all;
                            text_pos = [0.5, 0.75];
                            text_color = 'm';
                            mdl = fitlm(X, Y);
                            y_fit = predict(mdl, x_fit');
                            plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                            if plot_stat
                                if bootstrap
                                    stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                    r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                    p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                    text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                else
                                    [r, p] = corr(X', Y', 'Type', 'Pearson');
                                    text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                end
                            end
    
                            % text(0.5, 0.85, ['k = ', num2str(mdl_all.Coefficients.Estimate(2))], 'Units', 'normalized');
                            % text(0.5, 0.8, ['b = ', num2str(mdl_all.Coefficients.Estimate(1))], 'Units', 'normalized');
                        end

                        xlabel(['Pre, within, ', data_type_labels{data_type}]);
                        ylabel(['Post, within, ', data_type_labels{data_type}]);
                        axis equal;

                        % origin axis
                        plot([0, 0], ylim, 'k--');
                        plot(xlim, [0, 0], 'k--');

                        if plot_legend
                            legend({"ACC", "VLPFC"}, "Location", "northwest");
                        end
                    end

                    %% 2. cortical neurons: pre across vs post across
                    plot_types = [1, 2, 4]; % 1: mean, 2: significant count, 3: significant mean, 4: significant fixed mean
                    for plot_idx = 1:3
                        nexttile(plot_idx+3);
                        data_type = plot_types(plot_idx);
                        hold on;
                        
                        % ACC
                        filter = cellfun(@(x) x.area == 1, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_ACC = cellfun(@(x) x.all_data(1, 3, 3, data_type), filtered_profile);
                        y_ACC = cellfun(@(x) x.all_data(2, 3, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_ACC) & ~isnan(y_ACC);
                        x_ACC = x_ACC(nan_filter);
                        y_ACC = y_ACC(nan_filter);
                        scatter(x_ACC, y_ACC, 10, 'r', 'x', 'Tag', 'ACC');

                        % VLPFC
                        filter = cellfun(@(x) x.area == 3, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_VLPFC = cellfun(@(x) x.all_data(1, 1, 3, data_type), filtered_profile);
                        y_VLPFC = cellfun(@(x) x.all_data(2, 1, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_VLPFC) & ~isnan(y_VLPFC);
                        x_VLPFC = x_VLPFC(nan_filter);
                        y_VLPFC = y_VLPFC(nan_filter);
                        scatter(x_VLPFC, y_VLPFC, 10, 'b', '+', 'Tag', 'VLPFC');

                        x_all = [x_ACC, x_VLPFC];
                        y_all = [y_ACC, y_VLPFC];

                        % linear regression
                        if max(x_all)>min(x_all)
                            x_fit = linspace(min(x_all), max(x_all), 100);
                            plot(x_fit, x_fit, 'k:', 'LineWidth', 1); % identity line
    
                            if length(x_ACC) > 1
                                X = x_ACC;
                                Y = y_ACC;
                                text_pos = [0.5, 0.95];
                                text_color = 'r';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            if length(x_VLPFC) > 1
                                X = x_VLPFC;
                                Y = y_VLPFC;
                                text_pos = [0.5, 0.85];
                                text_color = 'b';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            X = x_all;
                            Y = y_all;
                            text_pos = [0.5, 0.75];
                            text_color = 'm';
                            mdl = fitlm(X, Y);
                            y_fit = predict(mdl, x_fit');
                            plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                            if plot_stat
                                if bootstrap
                                    stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                    r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                    p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                    text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                else
                                    [r, p] = corr(X', Y', 'Type', 'Pearson');
                                    text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                end
                            end
    
                            % text(0.5, 0.85, ['k = ', num2str(mdl_all.Coefficients.Estimate(2))], 'Units', 'normalized');
                            % text(0.5, 0.8, ['b = ', num2str(mdl_all.Coefficients.Estimate(1))], 'Units', 'normalized');
                        end

                        xlabel(['Pre, across, ', data_type_labels{data_type}]);
                        ylabel(['Post, across, ', data_type_labels{data_type}]);
                        axis equal;

                        % origin axis
                        plot([0, 0], ylim, 'k--');
                        plot(xlim, [0, 0], 'k--');

                        if plot_legend
                            legend({"ACC", "VLPFC"}, "Location", "northwest");
                        end
                    end

                    %% 3. cortical neurons: pre connection with thalamus vs within difference
                    plot_types = [1, 2, 4]; % 1: mean, 2: significant count, 3: significant mean, 4: significant fixed mean
                    for plot_idx = 1:3
                        nexttile(plot_idx+6);
                        data_type = plot_types(plot_idx);
                        hold on;

                        % ACC
                        filter = cellfun(@(x) x.area == 1, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_ACC = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                        y_ACC = cellfun(@(x) x.all_data_diff(1, 1, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_ACC) & ~isnan(y_ACC);
                        x_ACC = x_ACC(nan_filter);
                        y_ACC = y_ACC(nan_filter);
                        scatter(x_ACC, y_ACC, 10, 'r', 'x', 'Tag', 'ACC');

                        % VLPFC
                        filter = cellfun(@(x) x.area == 3, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_VLPFC = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                        y_VLPFC = cellfun(@(x) x.all_data_diff(1, 3, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_VLPFC) & ~isnan(y_VLPFC);
                        x_VLPFC = x_VLPFC(nan_filter);
                        y_VLPFC = y_VLPFC(nan_filter);
                        scatter(x_VLPFC, y_VLPFC, 10, 'b', '+', 'Tag', 'VLPFC');

                        x_all = [x_ACC, x_VLPFC];
                        y_all = [y_ACC, y_VLPFC];

                        % linear regression
                        if max(x_all)>min(x_all)
                            x_fit = linspace(min(x_all), max(x_all), 100);
                            plot(x_fit, x_fit, 'k:', 'LineWidth', 1); % identity line
    
                            if length(x_ACC) > 1
                                X = x_ACC;
                                Y = y_ACC;
                                text_pos = [0.5, 0.95];
                                text_color = 'r';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            if length(x_VLPFC) > 1
                                X = x_VLPFC;
                                Y = y_VLPFC;
                                text_pos = [0.5, 0.85];
                                text_color = 'b';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            X = x_all;
                            Y = y_all;
                            text_pos = [0.5, 0.75];
                            text_color = 'm';
                            mdl = fitlm(X, Y);
                            y_fit = predict(mdl, x_fit');
                            plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                            if plot_stat
                                if bootstrap
                                    stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                    r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                    p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                    text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                else
                                    [r, p] = corr(X', Y', 'Type', 'Pearson');
                                    text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                end
                            end
    
                            % text(0.5, 0.85, ['k = ', num2str(mdl_all.Coefficients.Estimate(2))], 'Units', 'normalized');
                            % text(0.5, 0.8, ['b = ', num2str(mdl_all.Coefficients.Estimate(1))], 'Units', 'normalized');
                        end

                        xlabel(['Pre, J to Th, ', data_type_labels{data_type}]);
                        ylabel(['Post - Pre, within, ', data_type_labels{data_type}]);
                        axis equal;

                        % origin axis
                        plot([0, 0], ylim, 'k--');
                        plot(xlim, [0, 0], 'k--');

                        if plot_legend
                            legend({"ACC", "VLPFC"}, "Location", "northwest");
                        end
                    end

                    %% 4. cortical neurons: pre connection with thalamus vs across difference
                    plot_types = [1, 2, 4]; % 1: mean, 2: significant count, 3: significant mean, 4: significant fixed mean
                    for plot_idx = 1:3
                        nexttile(plot_idx+9);
                        data_type = plot_types(plot_idx);
                        hold on;
                        
                        % ACC
                        filter = cellfun(@(x) x.area == 1, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_ACC = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                        y_ACC = cellfun(@(x) x.all_data_diff(1, 3, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_ACC) & ~isnan(y_ACC);
                        x_ACC = x_ACC(nan_filter);
                        y_ACC = y_ACC(nan_filter);
                        scatter(x_ACC, y_ACC, 10, 'r', 'x', 'Tag', 'ACC');

                        % VLPFC
                        filter = cellfun(@(x) x.area == 3, cell_profile);
                        filtered_profile = cell_profile(filter);
                        x_VLPFC = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                        y_VLPFC = cellfun(@(x) x.all_data_diff(1, 1, 3, data_type), filtered_profile);
                        nan_filter = ~isnan(x_VLPFC) & ~isnan(y_VLPFC);
                        x_VLPFC = x_VLPFC(nan_filter);
                        y_VLPFC = y_VLPFC(nan_filter);
                        scatter(x_VLPFC, y_VLPFC, 10, 'b', '+', 'Tag', 'VLPFC');

                        x_all = [x_ACC, x_VLPFC];
                        y_all = [y_ACC, y_VLPFC];

                        % linear regression
                        if max(x_all)>min(x_all)
                            x_fit = linspace(min(x_all), max(x_all), 100);
                            plot(x_fit, x_fit, 'k:', 'LineWidth', 1); % identity line
    
                            if length(x_ACC) > 1
                                X = x_ACC;
                                Y = y_ACC;
                                text_pos = [0.5, 0.95];
                                text_color = 'r';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            if length(x_VLPFC) > 1
                                X = x_VLPFC;
                                Y = y_VLPFC;
                                text_pos = [0.5, 0.85];
                                text_color = 'b';
                                mdl = fitlm(X, Y);
                                y_fit = predict(mdl, x_fit');
                                plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                                if plot_stat
                                    if bootstrap
                                        stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                        r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                        p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                        text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                        text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    else
                                        [r, p] = corr(X', Y', 'Type', 'Pearson');
                                        text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                    end
                                end
                            end
    
                            X = x_all;
                            Y = y_all;
                            text_pos = [0.5, 0.75];
                            text_color = 'm';
                            mdl = fitlm(X, Y);
                            y_fit = predict(mdl, x_fit');
                            plot(x_fit, y_fit, [text_color, '--'], 'LineWidth', 1);
                            if plot_stat
                                if bootstrap
                                    stat_boot = bootstrp(1000, @corr_merged, X', Y');
                                    r_CI = prctile(stat_boot(:, 1), [2.5, 97.5]);
                                    p_CI = prctile(stat_boot(:, 2), [2.5, 97.5]);
                                    text(text_pos(1), text_pos(2), ['rho=[', num2str(round(r_CI(1), 3)), ', ', num2str(round(r_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                    text(text_pos(1), text_pos(2)-0.05, ['p=[', num2str(round(p_CI(1), 3)), ', ', num2str(round(p_CI(2), 3)), ']'], 'Units', 'normalized', 'Color', text_color);
                                else
                                    [r, p] = corr(X', Y', 'Type', 'Pearson');
                                    text(text_pos(1), text_pos(2), ['rho=', num2str(round(r, 3)), ', p=', num2str(round(p, 3))], 'Units', 'normalized', 'Color', text_color);
                                end
                            end
                            % text(0.5, 0.85, ['k = ', num2str(mdl_all.Coefficients.Estimate(2))], 'Units', 'normalized');
                            % text(0.5, 0.8, ['b = ', num2str(mdl_all.Coefficients.Estimate(1))], 'Units', 'normalized');
                        end

                        xlabel(['Pre, J to Th, ', data_type_labels{data_type}]);
                        ylabel(['Post - Pre, across, ', data_type_labels{data_type}]);
                        axis equal;

                        % origin axis
                        plot([0, 0], ylim, 'k--');
                        plot(xlim, [0, 0], 'k--');

                        if plot_legend
                            legend({"ACC", "VLPFC"}, "Location", "northwest");
                        end
                    end


                    %% save
                    if session_idx == 0
                        title(t, [session_type, ' All sessions, ', states{state_idx}, ', ', aligns{align_idx}, ', Kernel ', num2str(kernel_idx), ', ', sign_filter]);
                    else
                        title(t, [session_type, 'Session ', num2str(session_idx), ', ', states{state_idx}, ', ', aligns{align_idx}, ', Kernel ', num2str(kernel_idx), ', ', sign_filter]);
                    end

                    save_folder = [root_path, 'figures/network_io/', session_type, '/', sign_filter, '/diff'];
                    if ~exist(save_folder, 'dir')
                        mkdir(save_folder);
                    end
                    save_name = [save_folder, '/network_io_', states{state_idx}, '_', aligns{align_idx}, '_', num2str(session_idx), '_', num2str(kernel_idx), '.png'];
                    if bootstrap && plot_stat
                        save_name = [save_name(1:end-4), '_bootstrap', save_name(end-3:end)];
                    end
                    if ~plot_stat
                        save_name = [save_name(1:end-4), '_no_stat', save_name(end-3:end)];
                    end
                    % saveas(f, save_name);
                    print(f, save_name, '-dpng');
                    
                    close(f);

                end
                fprintf("Done.\n");
            end
        end
    end
end

function stat=corr_merged(x, y)
    [r, p] = corr(x, y, 'Type', 'Pearson');
    stat = [r, p];
end
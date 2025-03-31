% Use last aligned time bins

% session_type = 'Muscimol';
session_type = 'Saline';
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
sign_filter = 'raw'; % 'abs', 'pos', 'neg', or 'raw'

if strcmp(session_type, 'Muscimol')
    n_session = 10;
else
    n_session = 5;
end

n_conn_kernel = 3;
J_data = cell(3, 3, n_conn_kernel, n_states, 1, n_session, 2); % (area i, area j, kernel, state, align, session, prepost), each cell is a n_area_i x n_area_j matrix
J_data_err = cell(3, 3, n_conn_kernel, n_states, 1, n_session, 2);
area_size = zeros(3, 3, n_session);

all_J = [];
J_state = [];
J_session = [];
J_kernel = [];
J_i = [];
J_j = [];

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
                    cell_profile{i}.sig_count = zeros(2, 3, 3); % number of significant J to an area
                    cell_profile{i}.sig_mean = zeros(2, 3, 3); % mean of significant J to an area
                    cell_profile{i}.all_data = zeros(2, 3, 3, 3); % all J data

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

                        end
                    end
                    cell_profile{i}.all_data(:, :, :, 1) = cell_profile{i}.mean(:, :, :);
                    cell_profile{i}.all_data(:, :, :, 2) = cell_profile{i}.sig_count(:, :, :);
                    cell_profile{i}.all_data(:, :, :, 3) = cell_profile{i}.sig_mean(:, :, :);

                    % for cortical neurons, calculate post - pre for within and across area connections
                    if cell_profile{i}.area == 1 || cell_profile{i}.area == 3
                        cell_profile{i}.mean_diff = cell_profile{i}.mean(2, :, :) - cell_profile{i}.mean(1, :, :);
                        cell_profile{i}.sig_count_diff = cell_profile{i}.sig_count(2, :, :) - cell_profile{i}.sig_count(1, :, :);
                        cell_profile{i}.sig_mean_diff = cell_profile{i}.sig_mean(2, :, :) - cell_profile{i}.sig_mean(1, :, :);
                        cell_profile{i}.all_data_diff = cell_profile{i}.all_data(2, :, :, :) - cell_profile{i}.all_data(1, :, :, :);
                    end
                end
                all_profile = [all_profile, cell_profile];
            end

            %% plot
            for session_idx = 0:n_session
                %% Figure 1: Pre
                f = figure("Position", [0, 0, 1200, 1200], "Visible", "off");
                t = tiledlayout(3, 3, "TileSpacing", "compact", "Padding", "compact");

                if session_idx == 0
                    cell_profile = all_profile;
                else
                    filter = cellfun(@(x) x.session == session_idx, all_profile);
                    cell_profile = all_profile(filter);
                end

                data_type_labels = {'Mean J', 'Significant J Count', 'Mean of Significant J'};

                %% 1. cortical neurons: pre connection with thalamus vs pre within area connection
                for data_type = 1:3
                    nexttile(data_type);
                    hold on;
                    x_all = [];
                    y_all = [];
                    % ACC
                    filter = cellfun(@(x) x.area == 1, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data(1, 1, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'r', 'filled', 'Tag', 'ACC');

                    % VLPFC
                    filter = cellfun(@(x) x.area == 3, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data(1, 3, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'b', 'filled', 'Tag', 'VLPFC');

                    % linear regression
                    mdl = fitlm(x_all, y_all);
                    x_fit = linspace(min(x_all), max(x_all), 100);
                    y_fit = predict(mdl, x_fit');
                    plot(x_fit, y_fit, 'm--', 'LineWidth', 1);
                    text(0.5, 0.9, ['R^2 = ', num2str(mdl.Rsquared.Ordinary)], 'Units', 'normalized');

                    xlabel(['Pre, connection with Thalamus, ', data_type_labels{data_type}]);
                    ylabel(['Pre, within area connection, ', data_type_labels{data_type}]);
                    axis equal;

                    % origin axis
                    plot([0, 0], ylim, 'k--');
                    plot(xlim, [0, 0], 'k--');
                end

                %% 2. cortical neurons: pre connection with thalamus vs pre across area connection
                for data_type = 1:3
                    nexttile(data_type+3);
                    hold on;
                    x_all = [];
                    y_all = [];
                    % ACC
                    filter = cellfun(@(x) x.area == 1, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data(1, 3, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'r', 'filled', 'Tag', 'ACC');

                    % VLPFC
                    filter = cellfun(@(x) x.area == 3, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data(1, 1, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'b', 'filled', 'Tag', 'VLPFC');

                    % linear regression
                    mdl = fitlm(x_all, y_all);
                    x_fit = linspace(min(x_all), max(x_all), 100);
                    y_fit = predict(mdl, x_fit');
                    plot(x_fit, y_fit, 'm--', 'LineWidth', 1);
                    text(0.5, 0.9, ['R^2 = ', num2str(mdl.Rsquared.Ordinary)], 'Units', 'normalized');

                    xlabel(['Pre, connection with Thalamus, ', data_type_labels{data_type}]);
                    ylabel(['Pre, across area connection, ', data_type_labels{data_type}]);
                    axis equal;

                    % origin axis
                    plot([0, 0], ylim, 'k--');
                    plot(xlim, [0, 0], 'k--');
                end

                %% 3. cortical neurons: pre within area connection vs pre across area connection
                for data_type = 1:3
                    nexttile(data_type+6);
                    hold on;
                    x_all = [];
                    y_all = [];
                    % ACC
                    filter = cellfun(@(x) x.area == 1, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(1, 1, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data(1, 3, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'r', 'filled', 'Tag', 'ACC');

                    % VLPFC
                    filter = cellfun(@(x) x.area == 3, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(1, 3, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data(1, 1, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'b', 'filled', 'Tag', 'VLPFC');

                    % linear regression
                    mdl = fitlm(x_all, y_all);
                    x_fit = linspace(min(x_all), max(x_all), 100);
                    y_fit = predict(mdl, x_fit');
                    plot(x_fit, y_fit, 'm--', 'LineWidth', 1);
                    text(0.5, 0.9, ['R^2 = ', num2str(mdl.Rsquared.Ordinary)], 'Units', 'normalized');

                    xlabel(['Pre, within area connection, ', data_type_labels{data_type}]);
                    ylabel(['Pre, across area connection, ', data_type_labels{data_type}]);
                    axis equal;

                    % origin axis
                    plot([0, 0], ylim, 'k--');
                    plot(xlim, [0, 0], 'k--');
                end


                %% save
                if session_idx == 0
                    title(t, ['All sessions, ', states{state_idx}, ', ', aligns{align_idx}, ', Kernel ', num2str(kernel_idx)]);
                else
                    title(t, ['Session ', num2str(session_idx), ', ', states{state_idx}, ', ', aligns{align_idx}, ', Kernel ', num2str(kernel_idx)]);
                end

                save_folder = [root_path, 'figures/network_io/raw/pre'];
                if ~exist(save_folder, 'dir')
                    mkdir(save_folder);
                end
                save_name = [save_folder, '/network_io_', states{state_idx}, '_', aligns{align_idx}, '_', num2str(session_idx), '_', num2str(kernel_idx), '.png'];
                saveas(f, save_name);
                close(f);

                %% Figure 2: Post
                f = figure("Position", [0, 0, 1200, 1200], "Visible", "off");
                t = tiledlayout(3, 3, "TileSpacing", "compact", "Padding", "compact");

                if session_idx == 0
                    cell_profile = all_profile;
                else
                    filter = cellfun(@(x) x.session == session_idx, all_profile);
                    cell_profile = all_profile(filter);
                end

                data_type_labels = {'Mean J', 'Significant J Count', 'Mean of Significant J'};

                %% 1. cortical neurons: pre connection with thalamus vs post within area connection
                for data_type = 1:3
                    nexttile(data_type);
                    hold on;
                    x_all = [];
                    y_all = [];
                    % ACC
                    filter = cellfun(@(x) x.area == 1, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data(2, 1, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'r', 'filled', 'Tag', 'ACC');

                    % VLPFC
                    filter = cellfun(@(x) x.area == 3, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data(2, 3, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'b', 'filled', 'Tag', 'VLPFC');

                    % linear regression
                    mdl = fitlm(x_all, y_all);
                    x_fit = linspace(min(x_all), max(x_all), 100);
                    y_fit = predict(mdl, x_fit');
                    plot(x_fit, y_fit, 'm--', 'LineWidth', 1);
                    text(0.5, 0.9, ['R^2 = ', num2str(mdl.Rsquared.Ordinary)], 'Units', 'normalized');

                    xlabel(['Pre, connection with Thalamus, ', data_type_labels{data_type}]);
                    ylabel(['Post, within area connection, ', data_type_labels{data_type}]);
                    axis equal;

                    % origin axis
                    plot([0, 0], ylim, 'k--');
                    plot(xlim, [0, 0], 'k--');
                end

                %% 2. cortical neurons: pre connection with thalamus vs post across area connection
                for data_type = 1:3
                    nexttile(data_type+3);
                    hold on;
                    x_all = [];
                    y_all = [];
                    % ACC
                    filter = cellfun(@(x) x.area == 1, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data(2, 3, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'r', 'filled', 'Tag', 'ACC');

                    % VLPFC
                    filter = cellfun(@(x) x.area == 3, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data(2, 1, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'b', 'filled', 'Tag', 'VLPFC');

                    % linear regression
                    mdl = fitlm(x_all, y_all);
                    x_fit = linspace(min(x_all), max(x_all), 100);
                    y_fit = predict(mdl, x_fit');
                    plot(x_fit, y_fit, 'm--', 'LineWidth', 1);
                    text(0.5, 0.9, ['R^2 = ', num2str(mdl.Rsquared.Ordinary)], 'Units', 'normalized');

                    xlabel(['Pre, connection with Thalamus, ', data_type_labels{data_type}]);
                    ylabel(['Post, across area connection, ', data_type_labels{data_type}]);
                    axis equal;

                    % origin axis
                    plot([0, 0], ylim, 'k--');
                    plot(xlim, [0, 0], 'k--');
                end

                %% 3. cortical neurons: post within area connection vs post across area connection
                for data_type = 1:3
                    nexttile(data_type+6);
                    hold on;
                    x_all = [];
                    y_all = [];
                    % ACC
                    filter = cellfun(@(x) x.area == 1, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(2, 1, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data(2, 3, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'r', 'filled', 'Tag', 'ACC');

                    % VLPFC
                    filter = cellfun(@(x) x.area == 3, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(2, 3, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data(2, 1, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'b', 'filled', 'Tag', 'VLPFC');

                    % linear regression
                    mdl = fitlm(x_all, y_all);
                    x_fit = linspace(min(x_all), max(x_all), 100);
                    y_fit = predict(mdl, x_fit');
                    plot(x_fit, y_fit, 'm--', 'LineWidth', 1);
                    text(0.5, 0.9, ['R^2 = ', num2str(mdl.Rsquared.Ordinary)], 'Units', 'normalized');

                    xlabel(['Pre, within area connection, ', data_type_labels{data_type}]);
                    ylabel(['Post, across area connection, ', data_type_labels{data_type}]);
                    axis equal;

                    % origin axis
                    plot([0, 0], ylim, 'k--');
                    plot(xlim, [0, 0], 'k--');
                end


                %% save
                if session_idx == 0
                    title(t, ['All sessions, ', states{state_idx}, ', ', aligns{align_idx}, ', Kernel ', num2str(kernel_idx)]);
                else
                    title(t, ['Session ', num2str(session_idx), ', ', states{state_idx}, ', ', aligns{align_idx}, ', Kernel ', num2str(kernel_idx)]);
                end

                save_folder = [root_path, 'figures/network_io/raw/post'];
                if ~exist(save_folder, 'dir')
                    mkdir(save_folder);
                end
                save_name = [save_folder, '/network_io_', states{state_idx}, '_', aligns{align_idx}, '_', num2str(session_idx), '_', num2str(kernel_idx), '.png'];
                saveas(f, save_name);
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

                data_type_labels = {'Mean J', 'Significant J Count', 'Mean of Significant J'};

                %% 1. cortical neurons: pre within area connection vs post within area connection
                for data_type = 1:3
                    nexttile(data_type);
                    hold on;
                    x_all = [];
                    y_all = [];
                    % ACC
                    filter = cellfun(@(x) x.area == 1, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(1, 1, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data(2, 1, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'r', 'filled', 'Tag', 'ACC');

                    % VLPFC
                    filter = cellfun(@(x) x.area == 3, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(1, 3, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data(2, 3, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'b', 'filled', 'Tag', 'VLPFC');

                    % linear regression
                    mdl = fitlm(x_all, y_all);
                    x_fit = linspace(min(x_all), max(x_all), 100);
                    y_fit = predict(mdl, x_fit');
                    plot(x_fit, y_fit, 'm--', 'LineWidth', 1);
                    text(0.5, 0.9, ['R^2 = ', num2str(mdl.Rsquared.Ordinary)], 'Units', 'normalized');

                    xlabel(['Pre, within area connection, ', data_type_labels{data_type}]);
                    ylabel(['Post, within area connection, ', data_type_labels{data_type}]);
                    axis equal;

                    % origin axis
                    plot([0, 0], ylim, 'k--');
                    plot(xlim, [0, 0], 'k--');
                end

                %% 2. cortical neurons: pre across area connection vs post across area connection
                for data_type = 1:3
                    nexttile(data_type+3);
                    hold on;
                    x_all = [];
                    y_all = [];
                    % ACC
                    filter = cellfun(@(x) x.area == 1, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(1, 3, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data(2, 3, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'r', 'filled', 'Tag', 'ACC');

                    % VLPFC
                    filter = cellfun(@(x) x.area == 3, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(1, 1, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data(2, 1, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'b', 'filled', 'Tag', 'VLPFC');

                    % linear regression
                    mdl = fitlm(x_all, y_all);
                    x_fit = linspace(min(x_all), max(x_all), 100);
                    y_fit = predict(mdl, x_fit');
                    plot(x_fit, y_fit, 'm--', 'LineWidth', 1);
                    text(0.5, 0.9, ['R^2 = ', num2str(mdl.Rsquared.Ordinary)], 'Units', 'normalized');

                    xlabel(['Pre, across area connection, ', data_type_labels{data_type}]);
                    ylabel(['Post, across area connection, ', data_type_labels{data_type}]);
                    axis equal;

                    % origin axis
                    plot([0, 0], ylim, 'k--');
                    plot(xlim, [0, 0], 'k--');
                end

                %% 3. cortical neurons: pre connection with thalamus vs within area connection difference
                for data_type = 1:3
                    nexttile(data_type+6);
                    hold on;
                    x_all = [];
                    y_all = [];
                    % ACC
                    filter = cellfun(@(x) x.area == 1, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data_diff(1, 1, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'r', 'filled', 'Tag', 'ACC');

                    % VLPFC
                    filter = cellfun(@(x) x.area == 3, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data_diff(1, 3, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'b', 'filled', 'Tag', 'VLPFC');

                    % linear regression
                    mdl = fitlm(x_all, y_all);
                    x_fit = linspace(min(x_all), max(x_all), 100);
                    y_fit = predict(mdl, x_fit');
                    plot(x_fit, y_fit, 'm--', 'LineWidth', 1);
                    text(0.5, 0.9, ['R^2 = ', num2str(mdl.Rsquared.Ordinary)], 'Units', 'normalized');

                    xlabel(['Pre, connection with Thalamus, ', data_type_labels{data_type}]);
                    ylabel(['Post - Pre, within area connection, ', data_type_labels{data_type}]);
                    axis equal;

                    % origin axis
                    plot([0, 0], ylim, 'k--');
                    plot(xlim, [0, 0], 'k--');
                end

                %% 4. cortical neurons: pre connection with thalamus vs across area connection difference
                for data_type = 1:3
                    nexttile(data_type+9);
                    hold on;
                    x_all = [];
                    y_all = [];
                    % ACC
                    filter = cellfun(@(x) x.area == 1, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data_diff(1, 3, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'r', 'filled', 'Tag', 'ACC');

                    % VLPFC
                    filter = cellfun(@(x) x.area == 3, cell_profile);
                    filtered_profile = cell_profile(filter);
                    x = cellfun(@(x) x.all_data(1, 2, 3, data_type), filtered_profile);
                    y = cellfun(@(x) x.all_data_diff(1, 1, 3, data_type), filtered_profile);
                    x_all = [x_all, x];
                    y_all = [y_all, y];
                    scatter(x, y, 10, 'b', 'filled', 'Tag', 'VLPFC');

                    % linear regression
                    mdl = fitlm(x_all, y_all);
                    x_fit = linspace(min(x_all), max(x_all), 100);
                    y_fit = predict(mdl, x_fit');
                    plot(x_fit, y_fit, 'm--', 'LineWidth', 1);
                    text(0.5, 0.9, ['R^2 = ', num2str(mdl.Rsquared.Ordinary)], 'Units', 'normalized');

                    xlabel(['Pre, connection with Thalamus, ', data_type_labels{data_type}]);
                    ylabel(['Post - Pre, across area connection, ', data_type_labels{data_type}]);
                    axis equal;

                    % origin axis
                    plot([0, 0], ylim, 'k--');
                    plot(xlim, [0, 0], 'k--');
                end


                %% save
                if session_idx == 0
                    title(t, ['All sessions, ', states{state_idx}, ', ', aligns{align_idx}, ', Kernel ', num2str(kernel_idx)]);
                else
                    title(t, ['Session ', num2str(session_idx), ', ', states{state_idx}, ', ', aligns{align_idx}, ', Kernel ', num2str(kernel_idx)]);
                end

                save_folder = [root_path, 'figures/network_io/raw/diff'];
                if ~exist(save_folder, 'dir')
                    mkdir(save_folder);
                end
                save_name = [save_folder, '/network_io_', states{state_idx}, '_', aligns{align_idx}, '_', num2str(session_idx), '_', num2str(kernel_idx), '.png'];
                saveas(f, save_name);
                close(f);

            end

        end
    end
end
% Correlation: J and multiplied GLM beta.

% add function path
addpath("UtilFuncs");
addpath('UtilFuncs/TLS/');
addpath("UtilFuncs/HELPER_GENERAL");

root_path = '../';
kernel = 'DeltaPure';
reg = 'L2=2';
epoch = '2500';
area_names = {'ACC', 'Thalamus', 'VLPFC'};
area_type_names = {'ACC', 'VLPFC', 'Across area', 'Within area'};
% area_type_names = {'Across area', 'Within area'};
n_area_types = length(area_type_names);
filter_threshold = 1;

% % switches
% selected_tuning_label = 'Cho2Info_IxU';
% tuning_method = 'or'; % 'and' or 'or', for tuning filter
% beta_method = 'product'; % 'product' or 'sum', for beta product
% apply_tuning_filter = true;
% apply_J_filter = true; 
use_abs_beta = true;
use_abs_J = true; 

% build all possible combinations of switches
all_options = {};
regressors = {'IxU', 'ER', 'UNC', 'Sub'};
stages = {'BeforeCho', 'Cho2Info', 'Info2Rew'};
tuning_methods = {'and', 'or'};
beta_methods = {'product', 'sum'};
for reg_idx = 1:length(regressors)
    for stage_idx = 1:length(stages)
        for i = 1:16
            all_options{end+1}.tuning_label = [stages{stage_idx}, '_', regressors{reg_idx}];
            n = i;
            k = mod(n, 2);
            n = floor(n / 2);
            all_options{end}.tuning_method = tuning_methods{k+1};
            k = mod(n, 2);
            n = floor(n / 2);
            all_options{end}.beta_method = beta_methods{k+1};
            k = mod(n, 2);
            n = floor(n / 2);
            all_options{end}.apply_tuning_filter = logical(k);
            k = mod(n, 2);
            n = floor(n / 2);
            all_options{end}.apply_J_filter = logical(k);
            all_options{end}.option_idx = i;
        end
    end
end

n_options = length(all_options);
fprintf('Total options: %d\n', n_options);
total_time = 0;

for option_idx = 1:n_options
    tic;
    close all;
    % get switches
    selected_tuning_label = all_options{option_idx}.tuning_label;
    tuning_method = all_options{option_idx}.tuning_method;
    beta_method = all_options{option_idx}.beta_method;
    apply_tuning_filter = all_options{option_idx}.apply_tuning_filter;
    apply_J_filter = all_options{option_idx}.apply_J_filter;
    % use_abs_beta = all_options{option_idx}.use_abs_beta;
    % use_abs_J = all_options{option_idx}.use_abs_J;

    % combine switches to a short string
    switches = [selected_tuning_label, '_M'];
    if strcmp(tuning_method, 'and')
        switches = [switches, '&'];
    elseif strcmp(tuning_method, 'or')
        switches = [switches, '^'];
    end
    if strcmp(beta_method, 'product')
        switches = [switches, 'x'];
    elseif strcmp(beta_method, 'sum')
        switches = [switches, '+'];
    end
    switches = sprintf('%sF%02dA%02d', switches, apply_J_filter*10 + apply_tuning_filter, use_abs_beta*10 + use_abs_J);
    fprintf('Switches: %s\n', switches);

    %% load data
    % states = {'Task', 'RestOpen', 'RestClose', 'All'};
    states = {'Task', 'RestOpen', 'RestClose'};
    aligns = {'AlignLast'};
    prepost_names = {'Pre', 'Post'};
    session_types = {'Muscimol', 'Saline'};
    session_nums = [8, 5];
    used_sessions = {[1,4,5,6,7,8,9,10], [1,2,3,4,5]};
    n_states = length(states);
    state_pairs = {[1, 2], [1, 3], [2, 3]};
    n_pairs = length(state_pairs);
    max_n_session = max(session_nums);
    n_conn_kernel = 3;

    % (area i, area j, kernel, state, align, session, prepost, session_type), each cell is a n_area_i x n_area_j matrix
    J_data = cell(3, 3, n_conn_kernel, n_states, 1, max_n_session, 2, 2); 
    J_data_err = cell(3, 3, n_conn_kernel, n_states, 1, max_n_session, 2, 2);
    tuning = cell(3, 1, max_n_session, 2, 2); % (area, align, session, prepost, session_type)
    beta = cell(3, 1, max_n_session, 2, 2); % (area, align, session, prepost, session_type)

    area_size = zeros(3, 3, max_n_session);

    all_J = [];
    J_state = [];
    J_session = [];
    J_kernel = [];
    J_i = [];
    J_j = [];
    for session_type_idx = 1:2
        session_type = session_types{session_type_idx};
        prepost_all = {'Pre', 'Post'};
        n_session = session_nums(session_type_idx);
        used_session = used_sessions{session_type_idx};
        for prepost_idx = 1:2
            prepost = prepost_all{prepost_idx};
            for session_idx = used_session
                fprintf('Loading session %d\n', session_idx);

                for state_idx = 1:n_states
                    state = states{state_idx};
                    if strcmp(prepost, 'Pre')
                        % session_stage = [session_type, 'Pre', state, '_full'];
                        session_stage = [session_type, 'Pre', state, '_cortex'];
                    else
                        session_stage = [session_type, 'Post', state, '_cortex'];
                    end

                    for align_idx = 1:1
                        align = aligns{align_idx};
                        session_stage_full = [session_stage, '_', align];

                        % load model parameters
                        file_path = [root_path, 'GLM_model/', session_stage_full,...
                            '/GLM_', session_stage_full, '_', num2str(session_idx), '_',...
                            kernel, '_0_', reg, '_', epoch, '.mat'];
                        fprintf('Loading %s\n', session_stage_full);
                        load(file_path, "model_par", "n_PS_kernel", "kernel_len", "N", "model_err");
                        load([root_path, 'GLM_data/', session_stage_full,'/borders_', session_stage_full, '_', ...
                                num2str(session_idx),'.mat'], "borders");
                        
                        borders = [1, borders+0.5]; % area i is from borders(i) to borders(i+1)
                        n_area = length(borders) - 1;

                        % load tuning data (redundant for all states but never mind)
                        file_path = [root_path, 'GLM_data/', session_stage_full,...
                            '/tuning_', session_stage_full, '_', num2str(session_idx),'.mat'];
                        load(file_path, "taskID_labels", "tuning_session", "beta_session");

                        selected_tuning_label_full = [prepost, selected_tuning_label];
                        tuning_idx = find(strcmp(taskID_labels, selected_tuning_label_full), 1);
                        selected_tuning = tuning_session(:, tuning_idx);
                        selected_beta = beta_session(:, tuning_idx);

                        for i = 1:n_area
                            if n_area == 2 && i == 2
                                % if i == 2
                                    i_eff = 3;
                                else
                                    i_eff = i;
                            end
                            tuning{i_eff, 1, session_idx, prepost_idx, session_type_idx} = selected_tuning(borders(i):borders(i+1)-1);
                            beta{i_eff, 1, session_idx, prepost_idx, session_type_idx} = selected_beta(borders(i):borders(i+1)-1);
                        end
                
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

    %% calculate correlation
    for session_type_idx = 1:2
        session_type = session_types{session_type_idx};
        n_session = session_nums(session_type_idx);
        used_session = used_sessions{session_type_idx};
        for kernel_idx = 1:n_conn_kernel
            for state_idx = 1:n_states
                state = states{state_idx};

                f = figure('Position', [100, 100, 1200, n_area_types*400], 'Visible', 'off');
                t = tiledlayout(n_area_types, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
                title(t, [session_type, ', Kernel ', num2str(kernel_idx), ', ', state]);

                plot_data = zeros(2, 4, n_area_types); %(Mus/Sal, rho/k/b/DeltaJ, area type)
                bootstrap_result_all = cell(1, 4); % bootstrap results to get p and CI. Each cell = (bootstrap, Mus/Sal, rho/k/b/DeltaJ)

                for area_type_idx = 1:n_area_types  % ACC, VLPFC, across, within
                    % store data for bootstrapping
                    bootstrap_data = cell(1, 2); % pre and post
                    area_type = area_type_names{area_type_idx};
                    
                    for prepost_idx = 1:2
                        prepost = prepost_names{prepost_idx};

                        nexttile;
                        % plot: X: beta product, Y: J
                        X = [];
                        Y = [];
                        sig_X = [];
                        sig_Y = [];

                        % collect data for plotting
                        for session_idx = used_session
                            switch area_type
                                case 'ACC'
                                    beta_vec = beta{1, 1, session_idx, prepost_idx, session_type_idx}; % (n, 1)
                                    if strcmp(beta_method, 'product')
                                        beta_product = beta_vec * beta_vec'; % (n, n)
                                    elseif strcmp(beta_method, 'sum')
                                        beta_product = beta_vec + beta_vec'; % (n, n)
                                    end
                                    tuning_vec = tuning{1, 1, session_idx, prepost_idx, session_type_idx}; % (n, 1)
                                    if strcmp(tuning_method, 'and')
                                        tuning_filter = tuning_vec & tuning_vec'; % (n, n)
                                    elseif strcmp(tuning_method, 'or')
                                        tuning_filter = tuning_vec | tuning_vec'; % (n, n)
                                    end
                                    J = J_data{1, 1, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx}; % (n, n)
                                    J_err = J_data_err{1, 1, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx}; % (n, n)
                                    mask = ~eye(size(beta_product)); % non-diagonal elements

                                    X = [X; beta_product(mask)];
                                    Y = [Y; J(mask)];
                                    sig_X = [sig_X; tuning_filter(mask)];
                                    sig_Y = [sig_Y; abs(beta_product(mask)) > filter_threshold * J_err(mask)];

                                case 'VLPFC'
                                    beta_vec = beta{3, 1, session_idx, prepost_idx, session_type_idx}; % (n, 1)
                                    if strcmp(beta_method, 'product')
                                        beta_product = beta_vec * beta_vec'; % (n, n)
                                    elseif strcmp(beta_method, 'sum')
                                        beta_product = beta_vec + beta_vec'; % (n, n)
                                    end
                                    tuning_vec = tuning{3, 1, session_idx, prepost_idx, session_type_idx}; % (n, 1)
                                    if strcmp(tuning_method, 'and')
                                        tuning_filter = tuning_vec & tuning_vec'; % (n, n)
                                    elseif strcmp(tuning_method, 'or')
                                        tuning_filter = tuning_vec | tuning_vec'; % (n, n)
                                    end
                                    J = J_data{3, 3, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx}; % (n, n)
                                    J_err = J_data_err{3, 3, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx}; % (n, n)
                                    mask = ~eye(size(beta_product)); % non-diagonal elements

                                    X = [X; beta_product(mask)];
                                    Y = [Y; J(mask)];
                                    sig_X = [sig_X; tuning_filter(mask)];
                                    sig_Y = [sig_Y; abs(beta_product(mask)) > filter_threshold * J_err(mask)];

                                case 'Across area'
                                    beta_vec1 = beta{1, 1, session_idx, prepost_idx, session_type_idx}; % (n, 1)
                                    beta_vec2 = beta{3, 1, session_idx, prepost_idx, session_type_idx}; % (n, 1)
                                    if strcmp(beta_method, 'product')
                                        beta_product = beta_vec1 * beta_vec2'; % (n1, n2)
                                    elseif strcmp(beta_method, 'sum')
                                        beta_product = beta_vec1 + beta_vec2'; % (n1, n2)
                                    end
                                    tuning_vec1 = tuning{1, 1, session_idx, prepost_idx, session_type_idx}; % (n, 1)
                                    tuning_vec2 = tuning{3, 1, session_idx, prepost_idx, session_type_idx}; % (n, 1)
                                    if strcmp(tuning_method, 'and')
                                        tuning_filter = tuning_vec1 & tuning_vec2'; % (n1, n2)
                                    elseif strcmp(tuning_method, 'or')
                                        tuning_filter = tuning_vec1 | tuning_vec2'; % (n1, n2)
                                    end
                                    J = J_data{1, 3, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx}; % (n, n)
                                    J_err = J_data_err{1, 3, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx}; % (n, n)
                                    assert(all(size(J) == size(beta_product)), 'J and beta product size mismatch');
                                    X = [X; beta_product(:)];
                                    Y = [Y; J(:)];
                                    sig_X = [sig_X; tuning_filter(:)];
                                    sig_Y = [sig_Y; abs(beta_product(:)) > filter_threshold * J_err(:)];
                                    if strcmp(beta_method, 'product')
                                        beta_product = beta_vec2 * beta_vec1'; % (n2, n1)
                                    elseif strcmp(beta_method, 'sum')
                                        beta_product = beta_vec2 + beta_vec1'; % (n2, n1)
                                    end
                                    if strcmp(tuning_method, 'and')
                                        tuning_filter = tuning_vec2 & tuning_vec1'; % (n2, n1)
                                    elseif strcmp(tuning_method, 'or')
                                        tuning_filter = tuning_vec2 | tuning_vec1'; % (n2, n1)
                                    end
                                    J = J_data{3, 1, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx}; % (n, n)
                                    J_err = J_data_err{3, 1, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx}; % (n, n)
                                    assert(all(size(J) == size(beta_product)), 'J and beta product size mismatch');
                                    X = [X; beta_product(:)];
                                    Y = [Y; J(:)];
                                    sig_X = [sig_X; tuning_filter(:)];
                                    sig_Y = [sig_Y; abs(beta_product(:)) > filter_threshold * J_err(:)];

                                case 'Within area'
                                    beta_vec = beta{1, 1, session_idx, prepost_idx, session_type_idx}; % (n, 1)
                                    if strcmp(beta_method, 'product')
                                        beta_product = beta_vec * beta_vec'; % (n, n)
                                    elseif strcmp(beta_method, 'sum')
                                        beta_product = beta_vec + beta_vec'; % (n, n)
                                    end
                                    tuning_vec = tuning{1, 1, session_idx, prepost_idx, session_type_idx}; % (n, 1)
                                    if strcmp(tuning_method, 'and')
                                        tuning_filter = tuning_vec & tuning_vec'; % (n, n)
                                    elseif strcmp(tuning_method, 'or')
                                        tuning_filter = tuning_vec | tuning_vec'; % (n, n)
                                    end
                                    J = J_data{1, 1, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx}; % (n, n)
                                    J_err = J_data_err{1, 1, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx}; % (n, n)
                                    mask = ~eye(size(beta_product)); % non-diagonal elements
                                    X = [X; beta_product(mask)];
                                    Y = [Y; J(mask)];
                                    sig_X = [sig_X; tuning_filter(mask)];
                                    sig_Y = [sig_Y; abs(beta_product(mask)) > filter_threshold * J_err(mask)];
                                    beta_vec = beta{3, 1, session_idx, prepost_idx, session_type_idx}; % (n, 1)
                                    if strcmp(beta_method, 'product')
                                        beta_product = beta_vec * beta_vec'; % (n, n)
                                    elseif strcmp(beta_method, 'sum')
                                        beta_product = beta_vec + beta_vec'; % (n, n)
                                    end
                                    tuning_vec = tuning{3, 1, session_idx, prepost_idx, session_type_idx}; % (n, 1)
                                    if strcmp(tuning_method, 'and')
                                        tuning_filter = tuning_vec & tuning_vec'; % (n, n)
                                    elseif strcmp(tuning_method, 'or')
                                        tuning_filter = tuning_vec | tuning_vec'; % (n, n)
                                    end
                                    J = J_data{3, 3, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx}; % (n, n)
                                    J_err = J_data_err{3, 3, kernel_idx, state_idx, 1, session_idx, prepost_idx, session_type_idx}; % (n, n)
                                    mask = ~eye(size(beta_product)); % non-diagonal elements
                                    X = [X; beta_product(mask)];
                                    Y = [Y; J(mask)];
                                    sig_X = [sig_X; tuning_filter(mask)];
                                    sig_Y = [sig_Y; abs(beta_product(mask)) > filter_threshold * J_err(mask)];

                            end
                        end

                        % NaN filter
                        filter = ~isnan(X) & ~isnan(Y);
                        if apply_tuning_filter
                            filter = filter & sig_X; % only keep significant beta product
                        end
                        if apply_J_filter
                            filter = filter & sig_Y; % only keep significant J
                        end
                        X = X(filter);
                        Y = Y(filter);
                        sig_X = sig_X(filter);
                        sig_Y = sig_Y(filter);

                        if use_abs_beta
                            X = abs(X); % use absolute beta product
                        end
                        if use_abs_J
                            Y = abs(Y); % use absolute J
                        end

                        % calculate axes limits, equal axes
                        x_lim = [min(X), max(X)];
                        y_lim = [min(Y), max(Y)];
                        % if x_lim(2)-x_lim(1) < y_lim(2)-y_lim(1)
                        %     lim_diff = (y_lim(2) - y_lim(1)) - (x_lim(2) - x_lim(1));
                        %     x_lim = [x_lim(1)-lim_diff/2, x_lim(2)+lim_diff/2];
                        % else
                        %     lim_diff = (x_lim(2) - x_lim(1)) - (y_lim(2) - y_lim(1));
                        %     y_lim = [y_lim(1)-lim_diff/2, y_lim(2)+lim_diff/2];
                        % end

                        % % x_lim = [-2, 6];
                        % % y_lim = [-2, 6];
                        % x_lim = [-1.5, 3.5];
                        % y_lim = [-1.5, 3.5];

                        % ----------------------------------plot
                        hold on;
                        % axis equal; 
                        % scatter(X, Y, 10, 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.1);
                        % scatter(X, Y, 10, '+', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

                        % 2D histogram
                        N_bins = 100;
                        xedges = linspace(x_lim(1), x_lim(2), N_bins);
                        yedges = linspace(y_lim(1), y_lim(2), N_bins);
                        [N, ~, ~] = histcounts2(X, Y, xedges, yedges);
                        
                        % plot 2D histogram, white to red colormap
                        cmap = brewermap(N_bins, 'Reds');
                        % Normalize the first color to white
                        cmap = (cmap ./ cmap(1, :))*0.99;

                        N_log1p = log1p(N); % log scale
                        imagesc(xedges, yedges, N_log1p'); % transpose N for correct orientation
                        set(gca, 'YDir', 'normal'); % flip y-axis
                        % show box border
                        set(gca, 'Box', 'on', 'XColor', 'k', 'YColor', 'k', 'FontSize', 8, 'LineWidth', 1);
                        ctick_labels = [0, 1, 10, 100, 1000, 10000]; % colorbar tick labels
                        cticks = log1p(ctick_labels); % log scale
                        colormap(cmap);
                        colorbar('Location', 'eastoutside', 'Ticks', cticks, 'TickLabels', ctick_labels);
                        clim(log1p([0, max(N(:))]));

                        % set(gca, 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'FontSize', 8, 'Box', 'off');
                        % set(gca, 'XTick', linspace(x_lim(1), x_lim(2), 5), 'YTick', linspace(y_lim(1), y_lim(2), 5));
                        % set(gca, 'XTickLabel', round(linspace(x_lim(1), x_lim(2), 5), 2), 'YTickLabel', round(linspace(y_lim(1), y_lim(2), 5), 2));
                        % set(gca, 'TickLength', [0.02, 0.02], 'LineWidth', 1, 'FontSize', 8, 'FontName', 'Arial');
                        % set(gca, 'Position', [0.1, 0.1, 0.8, 0.8]); % adjust position to fit colorbar

                        % axes
                        % x_lim = get(gca, 'xlim');
                        % y_lim = get(gca, 'ylim');
                        max_lim = [min(x_lim(1), y_lim(1)), max(x_lim(2), y_lim(2))];
                        % plot(max_lim, max_lim, 'k--', 'LineWidth', 0.7);
                        plot(x_lim, [0 0], 'k--', 'LineWidth', 0.7);
                        plot([0 0], y_lim, 'k--', 'LineWidth', 0.7);

                        % correlation line, total least square fit
                        [Err, P] = fit_2D_data(X, Y, 'no');
                        k = P(1); b = P(2);
                        y_pred = k * x_lim + b;
                        plot(x_lim, y_pred, 'm-', 'LineWidth', 0.7);
                        xlim(x_lim);
                        ylim(y_lim);
                        hold off;

                        % pearson correlation
                        [R, P] = corrcoef(X, Y);

                        % show stats
                        text(0.1, 0.9, ['n = ', num2str(length(X))], 'Units', 'normalized', 'FontSize', 8, 'Color', 'k');
                        text(0.1, 0.82, ['R = ', num2str(R(1, 2), '%.3f')], 'Units', 'normalized', 'FontSize', 8, 'Color', 'k');
                        text(0.1, 0.74, ['p = ', num2str(P(1, 2), '%.3f')], 'Units', 'normalized', 'FontSize', 8, 'Color', 'k');
                        text(0.1, 0.66, ['y = ', num2str(b, '%.3f'), '+', num2str(k, '%.3f'), 'x'], 'Units', 'normalized', 'FontSize', 8, 'Color', 'k');

                        % labels
                        title([prepost, ', ', area_type_names{area_type_idx}]);
                        xlabel('Beta product');
                        ylabel('J');
                        % xlim([-1, 1]);
                        % ylim([-1, 1]);

                        % save plot data
                        plot_data(prepost_idx, :, area_type_idx) = [R(1, 2), k, b, mean(abs(Y - X))];

                        % save bootstrap data
                        bootstrap_data{prepost_idx}.X = X;
                        bootstrap_data{prepost_idx}.Y = Y;
                        bootstrap_data{prepost_idx}.sig_X = sig_X;
                        bootstrap_data{prepost_idx}.sig_Y = sig_Y;
                    end

                    % bootstrap
                    Pre_X = bootstrap_data{1}.X;
                    Post_X = bootstrap_data{2}.X;
                    Pre_Y = bootstrap_data{1}.Y;
                    Post_Y = bootstrap_data{2}.Y;

                    n_bootstrap = 1000;
                    bootstrap_results = zeros(n_bootstrap, 2, 4); % (bootstrap, Pre/Post, rho/k/b/DeltaJ)
                
                    for i=1:n_bootstrap
                        Pre_idx = randi(length(Pre_X), length(Pre_X), 1);
                        Post_idx = randi(length(Post_X), length(Post_X), 1);
                        Pre_X_boot = Pre_X(Pre_idx);
                        Post_X_boot = Post_X(Post_idx);
                        Pre_Y_boot = Pre_Y(Pre_idx);
                        Post_Y_boot = Post_Y(Post_idx);

                        % pearson correlation
                        [Pre_R_boot, ~] = corrcoef(Pre_X_boot, Pre_Y_boot);
                        [Post_R_boot, ~] = corrcoef(Post_X_boot, Post_Y_boot);
                        % total least square fit
                        [Pre_Err, Pre_P] = fit_2D_data(Pre_X_boot, Pre_Y_boot, 'no');
                        [Post_Err, Post_P] = fit_2D_data(Post_X_boot, Post_Y_boot, 'no');
                        % delta J
                        Pre_DeltaJ = mean(abs(Pre_Y_boot - Pre_X_boot));
                        Post_DeltaJ = mean(abs(Post_Y_boot - Post_X_boot));

                        % save results
                        bootstrap_results(i, :, 1) = [Pre_R_boot(1, 2), Post_R_boot(1, 2)];
                        bootstrap_results(i, :, 2) = [Pre_P(1), Post_P(1)];
                        bootstrap_results(i, :, 3) = [Pre_P(2), Post_P(2)];
                        bootstrap_results(i, :, 4) = [Pre_DeltaJ, Post_DeltaJ];
                    end
                    bootstrap_result_all{area_type_idx} = bootstrap_results;

                    r_diff = bootstrap_results(:, 1, 1) - bootstrap_results(:, 2, 1);
                    p_left = nnz(r_diff < 0) / n_bootstrap;
                    p_right = nnz(r_diff > 0) / n_bootstrap;
                    p = min(p_left, p_right) * 2; % two-tailed p-value
                    ci = prctile(r_diff, [2.5, 97.5]);
                    ci1 = prctile(bootstrap_results(:, 1, 1), [2.5, 97.5]);
                    ci2 = prctile(bootstrap_results(:, 2, 1), [2.5, 97.5]);

                    % show results
                    nexttile;
                    set(gca, 'XColor', 'none', 'YColor', 'none', 'XTick', [], 'YTick', []);

                    text(0.5, 0.7, ['p = ', num2str(p, '%.3f')], 'Units', 'normalized', 'FontSize', 12, 'Color', 'k','HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                    text(0.5, 0.6, ['Pre rho CI = [', num2str(ci1(1), '%.3f'), ', ', num2str(ci1(2), '%.3f'), ']'], 'Units', 'normalized', 'FontSize', 12, 'Color', 'k','HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                    text(0.5, 0.5, ['Post rho CI = [', num2str(ci2(1), '%.3f'), ', ', num2str(ci2(2), '%.3f'), ']'], 'Units', 'normalized', 'FontSize', 12, 'Color', 'k','HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                    text(0.5, 0.4, ['rho diff CI = [', num2str(ci(1), '%.3f'), ', ', num2str(ci(2), '%.3f'), ']'], 'Units', 'normalized', 'FontSize', 12, 'Color', 'k','HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                    title('Bootstrap p-value');
                    
                end
                % save figure
                file_folder = [root_path, 'figures/GLM/bootstrap/', switches,'/',];
                if ~exist(file_folder, 'dir')
                    mkdir(file_folder);
                end
                file_name = [file_folder, 'TuningCorr_', session_type, '_', state, '_', num2str(kernel_idx), '.png'];
                saveas(f, file_name);
                close(f);
                
                % Plot bootstrap bar plots
                f = figure('Position', [100, 100, 1200, 1200], 'Visible', 'off');
                t = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
                title(t, [session_type, ' Kernel ', num2str(kernel_idx), ', ', state]);

                data_types = {'rho', 'k(slope)', 'b(offset)', 'J diff'};
                for data_type_idx = 1:4
                    nexttile;
                    title(data_types{data_type_idx});
                    hold on;
                    for area_type_idx = 1:n_area_types
                        bootstrap_results = bootstrap_result_all{area_type_idx};
                        Pre_data = bootstrap_results(:, 1, data_type_idx);
                        Post_data = bootstrap_results(:, 2, data_type_idx);
                        
                        bar(area_type_idx-0.2, plot_data(1, data_type_idx, area_type_idx), 0.4, 'FaceColor', [0.4, 0.4, 0.4], 'EdgeColor', 'k', 'LineWidth', 1);
                        bar(area_type_idx+0.2, plot_data(2, data_type_idx, area_type_idx), 0.4, 'FaceColor', [0.85, 0.85, 0.85], 'EdgeColor', 'k', 'LineWidth', 1);

                        % plot CI
                        CI_Pre = prctile(Pre_data, [2.5, 97.5]);
                        CI_Post = prctile(Post_data, [2.5, 97.5]);
                        CI_Pre_center = mean(CI_Pre);
                        CI_Post_center = mean(CI_Post);
                        CI_Pre_width = CI_Pre(2) - CI_Pre(1);
                        CI_Post_width = CI_Post(2) - CI_Post(1);
                        errorbar(area_type_idx-0.2, CI_Pre_center, CI_Pre_width/2, 'k', 'LineStyle', 'none', 'LineWidth', 1);
                        errorbar(area_type_idx+0.2, CI_Post_center, CI_Post_width/2, 'k', 'LineStyle', 'none', 'LineWidth', 1);

                        % significance
                        p_Pre_left = sum(Pre_data < 0) / length(Pre_data);
                        p_Pre_right = sum(Pre_data > 0) / length(Pre_data);
                        p_Pre = min(p_Pre_left, p_Pre_right) * 2;
                        p_Post_left = sum(Post_data < 0) / length(Post_data);
                        p_Post_right = sum(Post_data > 0) / length(Post_data);
                        p_Post = min(p_Post_left, p_Post_right) * 2;

                        text_x = area_type_idx - 0.2; text_y = CI_Pre(2);
                        if p_Pre < 0.001
                            text(text_x, text_y, '***', 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                        elseif p_Pre < 0.01
                            text(text_x, text_y, '**', 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                        elseif p_Pre < 0.05
                            text(text_x, text_y, '*', 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                        end

                        text_x = area_type_idx + 0.2; text_y = CI_Post(2);
                        if p_Post < 0.001
                            text(text_x, text_y, '***', 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                        elseif p_Post < 0.01
                            text(text_x, text_y, '**', 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                        elseif p_Post < 0.05
                            text(text_x, text_y, '*', 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                        end

                        p_diff_left = sum(Pre_data - Post_data < 0) / length(Pre_data);
                        p_diff_right = sum(Pre_data - Post_data > 0) / length(Pre_data);
                        p_diff = min(p_diff_left, p_diff_right) * 2;

                        text_x = area_type_idx; text_y = max(CI_Pre(2), CI_Post(2))*1.2;
                        if p_diff < 0.001
                            text(text_x, text_y, '***', 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
                            plot([text_x-0.2, text_x+0.2], [text_y,text_y], 'k-', 'LineWidth', 1);
                        elseif p_diff < 0.01
                            text(text_x, text_y, '**', 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
                            plot([text_x-0.2, text_x+0.2], [text_y,text_y], 'k-', 'LineWidth', 1);
                        elseif p_diff < 0.05
                            text(text_x, text_y, '*', 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
                            plot([text_x-0.2, text_x+0.2], [text_y,text_y], 'k-', 'LineWidth', 1);
                        end
                    end
                    hold off;
                    legend({'Pre', 'Post'}, 'Location', 'northwest', 'FontSize', 12, 'Box', 'off');
                    ylabel('Value');
                    set(gca, 'XTick', 1:n_area_types, 'XTickLabel', area_type_names, 'FontSize', 12);
                    xlim([0.5, n_area_types+0.5]);
                end

                % save figure
                file_name = [file_folder, 'TuningCorr_hist_', session_type, '_', state, '_', num2str(kernel_idx), '.png'];
                saveas(f, file_name);
                close(f);

            end
        end
    end

    t = toc;
    fprintf('Option time: %.2f seconds\n', t);
    total_time = total_time + t;
    fprintf('Total time: %.2f seconds\n', total_time);
    eta_time = total_time / option_idx * (length(session_types) - option_idx);
    fprintf('%d/%d sessions done. eta: %d:%d:%d\n', ...
        option_idx, n_options, floor(eta_time/3600), ...
        floor(mod(eta_time, 3600)/60), ...
        mod(eta_time, 60));
    fprintf('----------------------------------------\n');
end
% pre-post stability analysis

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

%% load data
states = {'Task', 'RestOpen', 'RestClose', 'All'};
aligns = {'AlignLast'};
session_types = {'Muscimol', 'Saline'};
session_nums = [8, 5];
n_states = length(states);
n_session = 10;
n_conn_kernel = 3;
used_sessions = {[1,4,5,6,7,8,9,10], [1,2,3,4,5]};

% (area i, area j, kernel, state, align, session, prepost, session_type), each cell is a n_area_i x n_area_j matrix
J_data = cell(3, 3, n_conn_kernel, n_states, 1, n_session, 2, 2); 
J_data_err = cell(3, 3, n_conn_kernel, n_states, 1, n_session, 2, 2);
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
    n_session = session_nums(session_type_idx);
    used_session = used_sessions{session_type_idx};
    for prepost_idx = 1:2
        prepost = prepost_all{prepost_idx};
        for session_idx = used_session
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

%% main plot

significant_count = zeros(3, 3, n_area_types, n_conn_kernel, n_states, 2); % (X:pos/no/neg, Y:pos/no/neg, area_type, kernel, state, session type)

for kernel_idx = 1:3
    for state_idx = 1:n_states
        state = states{state_idx};
        f = figure('Position', [100, 100, 1200, n_area_types*400], 'Visible', 'off');
        t = tiledlayout(n_area_types, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
        title(t, [state, ', Kernel ', num2str(kernel_idx)]);

        plot_data = zeros(2, 4, n_area_types); %(Mus/Sal, rho/k/b/DeltaJ, area type)
        bootstrap_result_all = cell(1, 4); % bootstrap results to get p and CI. Each cell = (bootstrap, Mus/Sal, rho/k/b/DeltaJ)

        for area_type_idx = 1:n_area_types  % ACC, VLPFC, across, within
            % store data for bootstrapping
            bootstrap_data = cell(1, 2);
            area_type = area_type_names{area_type_idx};
            
            for session_type_idx = 1:2
                n_session = session_nums(session_type_idx);
                used_session = used_sessions{session_type_idx};

                nexttile;
                % plot: pre-post J correlation
                X = [];
                Y = [];
                sig_X = [];
                sig_Y = [];

                % collect data for plotting
                for session_idx = used_session
                    switch area_type
                        case 'ACC'
                            J_pre = J_data{1, 1, kernel_idx, state_idx, 1, session_idx, 1, session_type_idx};
                            J_post = J_data{1, 1, kernel_idx, state_idx, 1, session_idx, 2, session_type_idx};
                            J_pre_err = J_data_err{1, 1, kernel_idx, state_idx, 1, session_idx, 1, session_type_idx};
                            J_post_err = J_data_err{1, 1, kernel_idx, state_idx, 1, session_idx, 2, session_type_idx};
                            X = [X; J_pre(:)];
                            Y = [Y; J_post(:)];
                            sig_X = [sig_X; abs(J_pre(:)) > filter_threshold * J_pre_err(:)];
                            sig_Y = [sig_Y; abs(J_post(:)) > filter_threshold * J_post_err(:)];

                        case 'VLPFC'
                            J_pre = J_data{3, 3, kernel_idx, state_idx, 1, session_idx, 1, session_type_idx};
                            J_post = J_data{3, 3, kernel_idx, state_idx, 1, session_idx, 2, session_type_idx};
                            J_pre_err = J_data_err{3, 3, kernel_idx, state_idx, 1, session_idx, 1, session_type_idx};
                            J_post_err = J_data_err{3, 3, kernel_idx, state_idx, 1, session_idx, 2, session_type_idx};
                            X = [X; J_pre(:)];
                            Y = [Y; J_post(:)];
                            sig_X = [sig_X; abs(J_pre(:)) > filter_threshold * J_pre_err(:)];
                            sig_Y = [sig_Y; abs(J_post(:)) > filter_threshold * J_post_err(:)];

                        case 'Across area'
                            J_pre = J_data{1, 3, kernel_idx, state_idx, 1, session_idx, 1, session_type_idx};
                            J_post = J_data{1, 3, kernel_idx, state_idx, 1, session_idx, 2, session_type_idx};
                            J_pre_err = J_data_err{1, 3, kernel_idx, state_idx, 1, session_idx, 1, session_type_idx};
                            J_post_err = J_data_err{1, 3, kernel_idx, state_idx, 1, session_idx, 2, session_type_idx};
                            X = [X; J_pre(:)];
                            Y = [Y; J_post(:)];
                            sig_X = [sig_X; abs(J_pre(:)) > filter_threshold * J_pre_err(:)];
                            sig_Y = [sig_Y; abs(J_post(:)) > filter_threshold * J_post_err(:)];
                            J_pre = J_data{3, 1, kernel_idx, state_idx, 1, session_idx, 1, session_type_idx};
                            J_post = J_data{3, 1, kernel_idx, state_idx, 1, session_idx, 2, session_type_idx};
                            J_pre_err = J_data_err{3, 1, kernel_idx, state_idx, 1, session_idx, 1, session_type_idx};
                            J_post_err = J_data_err{3, 1, kernel_idx, state_idx, 1, session_idx, 2, session_type_idx};
                            X = [X; J_pre(:)];
                            Y = [Y; J_post(:)];
                            sig_X = [sig_X; abs(J_pre(:)) > filter_threshold * J_pre_err(:)];
                            sig_Y = [sig_Y; abs(J_post(:)) > filter_threshold * J_post_err(:)];

                        case 'Within area'
                            J_pre = J_data{1, 1, kernel_idx, state_idx, 1, session_idx, 1, session_type_idx};
                            J_post = J_data{1, 1, kernel_idx, state_idx, 1, session_idx, 2, session_type_idx};
                            J_pre_err = J_data_err{1, 1, kernel_idx, state_idx, 1, session_idx, 1, session_type_idx};
                            J_post_err = J_data_err{1, 1, kernel_idx, state_idx, 1, session_idx, 2, session_type_idx};
                            X = [X; J_pre(:)];
                            Y = [Y; J_post(:)];
                            sig_X = [sig_X; abs(J_pre(:)) > filter_threshold * J_pre_err(:)];
                            sig_Y = [sig_Y; abs(J_post(:)) > filter_threshold * J_post_err(:)];
                            J_pre = J_data{3, 3, kernel_idx, state_idx, 1, session_idx, 1, session_type_idx};
                            J_post = J_data{3, 3, kernel_idx, state_idx, 1, session_idx, 2, session_type_idx};
                            J_pre_err = J_data_err{3, 3, kernel_idx, state_idx, 1, session_idx, 1, session_type_idx};
                            J_post_err = J_data_err{3, 3, kernel_idx, state_idx, 1, session_idx, 2, session_type_idx};
                            X = [X; J_pre(:)];
                            Y = [Y; J_post(:)];
                            sig_X = [sig_X; abs(J_pre(:)) > filter_threshold * J_pre_err(:)];
                            sig_Y = [sig_Y; abs(J_post(:)) > filter_threshold * J_post_err(:)];

                    end
                end

                % NaN filter
                filter = ~isnan(X) & ~isnan(Y);
                X = X(filter);
                Y = Y(filter);
                sig_X = sig_X(filter);
                sig_Y = sig_Y(filter);

                % % significant and sign filter
                % filter = sig_X | sig_Y;
                % X = X(filter);
                % Y = Y(filter);
                % sig_X = sig_X(filter);
                % sig_Y = sig_Y(filter);

                % store significant count data
                for i=1:3
                    for j=1:3
                        filter = ones(size(X));
                        if i == 1
                            filter = filter & sig_X & X > 0; % positive pre
                        elseif i == 2
                            filter = filter & ~sig_X; % non-significant pre
                        elseif i == 3
                            filter = filter & sig_X & X < 0; % negative pre
                        end
                        if j == 1
                            filter = filter & sig_Y & Y > 0; % positive post
                        elseif j == 2
                            filter = filter & ~sig_Y; % non-significant post
                        elseif j == 3
                            filter = filter & sig_Y & Y < 0; % negative post
                        end
                        significant_count(i, j, area_type_idx, kernel_idx, state_idx, session_type_idx) = sum(filter);
                    end
                end

                % calculate axes limits, equal axes
                x_lim = [min(X), max(X)];
                y_lim = [min(Y), max(Y)];
                if x_lim(2)-x_lim(1) < y_lim(2)-y_lim(1)
                    lim_diff = (y_lim(2) - y_lim(1)) - (x_lim(2) - x_lim(1));
                    x_lim = [x_lim(1)-lim_diff/2, x_lim(2)+lim_diff/2];
                else
                    lim_diff = (x_lim(2) - x_lim(1)) - (y_lim(2) - y_lim(1));
                    y_lim = [y_lim(1)-lim_diff/2, y_lim(2)+lim_diff/2];
                end

                % x_lim = [-2, 6];
                % y_lim = [-2, 6];
                x_lim = [-1.5, 3.5];
                y_lim = [-1.5, 3.5];

                % ---plot
                hold on;
                axis equal; 
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
                cmap = cmap ./ cmap(1, :);

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

                % axes and y=x line
                % x_lim = get(gca, 'xlim');
                % y_lim = get(gca, 'ylim');
                max_lim = [min(x_lim(1), y_lim(1)), max(x_lim(2), y_lim(2))];
                plot(max_lim, max_lim, 'k--', 'LineWidth', 0.7);
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
                title([session_types{session_type_idx}, ', ', area_type_names{area_type_idx}]);
                xlabel('Pre');
                ylabel('Post');
                % xlim([-1, 1]);
                % ylim([-1, 1]);

                % save plot data
                plot_data(session_type_idx, :, area_type_idx) = [R(1, 2), k, b, mean(abs(Y - X))];

                % save bootstrap data
                bootstrap_data{session_type_idx}.Pre = X;
                bootstrap_data{session_type_idx}.Post = Y;
                bootstrap_data{session_type_idx}.sig_Pre = sig_X;
                bootstrap_data{session_type_idx}.sig_Post = sig_Y;
            end
            % bootstrap
            Mus_X = bootstrap_data{1}.Pre;
            Sal_X = bootstrap_data{2}.Pre;
            Mus_Y = bootstrap_data{1}.Post;
            Sal_Y = bootstrap_data{2}.Post;

            n_bootstrap = 1000;
            bootstrap_results = zeros(n_bootstrap, 2, 4); % (bootstrap, Mus/Sal, rho/k/b/DeltaJ)
        
            for i=1:n_bootstrap
                Mus_idx = randi(length(Mus_X), length(Mus_X), 1);
                Sal_idx = randi(length(Sal_X), length(Sal_X), 1);
                Mus_X_boot = Mus_X(Mus_idx);
                Sal_X_boot = Sal_X(Sal_idx);
                Mus_Y_boot = Mus_Y(Mus_idx);
                Sal_Y_boot = Sal_Y(Sal_idx);

                % pearson correlation
                [Mus_R_boot, ~] = corrcoef(Mus_X_boot, Mus_Y_boot);
                [Sal_R_boot, ~] = corrcoef(Sal_X_boot, Sal_Y_boot);
                % total least square fit
                [Mus_Err, Mus_P] = fit_2D_data(Mus_X_boot, Mus_Y_boot, 'no');
                [Sal_Err, Sal_P] = fit_2D_data(Sal_X_boot, Sal_Y_boot, 'no');
                % delta J
                Mus_DeltaJ = mean(abs(Mus_Y_boot - Mus_X_boot));
                Sal_DeltaJ = mean(abs(Sal_Y_boot - Sal_X_boot));

                % save results
                bootstrap_results(i, :, 1) = [Mus_R_boot(1, 2), Sal_R_boot(1, 2)];
                bootstrap_results(i, :, 2) = [Mus_P(1), Sal_P(1)];
                bootstrap_results(i, :, 3) = [Mus_P(2), Sal_P(2)];
                bootstrap_results(i, :, 4) = [Mus_DeltaJ, Sal_DeltaJ];
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
            text(0.5, 0.6, ['Muscimol rho CI = [', num2str(ci1(1), '%.3f'), ', ', num2str(ci1(2), '%.3f'), ']'], 'Units', 'normalized', 'FontSize', 12, 'Color', 'k','HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            text(0.5, 0.5, ['Saline rho CI = [', num2str(ci2(1), '%.3f'), ', ', num2str(ci2(2), '%.3f'), ']'], 'Units', 'normalized', 'FontSize', 12, 'Color', 'k','HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            text(0.5, 0.4, ['rho diff CI = [', num2str(ci(1), '%.3f'), ', ', num2str(ci(2), '%.3f'), ']'], 'Units', 'normalized', 'FontSize', 12, 'Color', 'k','HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            title('Bootstrap p-value');
            
        end
        % save figure
        file_folder = [root_path, 'figures/GLM/bootstrap/'];
        if ~exist(file_folder, 'dir')
            mkdir(file_folder);
        end
        file_name = [file_folder, 'Stability_', state, '_', num2str(kernel_idx), '.png'];
        saveas(f, file_name);
        close(f);

        % Plot bootstrap bar plots
        f = figure('Position', [100, 100, 1200, 1200], 'Visible', 'off');
        t = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
        title(t, [state, ', Kernel ', num2str(kernel_idx)]);

        data_types = {'rho', 'k(slope)', 'b(offset)', 'J diff'};
        for data_type_idx = 1:4
            nexttile;
            title(data_types{data_type_idx});
            hold on;
            for area_type_idx = 1:n_area_types
                bootstrap_results = bootstrap_result_all{area_type_idx};
                Mus_data = bootstrap_results(:, 1, data_type_idx);
                Sal_data = bootstrap_results(:, 2, data_type_idx);
                
                bar(area_type_idx-0.2, plot_data(1, data_type_idx, area_type_idx), 0.4, 'FaceColor', [0.4, 0.4, 0.4], 'EdgeColor', 'k', 'LineWidth', 1);
                bar(area_type_idx+0.2, plot_data(2, data_type_idx, area_type_idx), 0.4, 'FaceColor', [0.85, 0.85, 0.85], 'EdgeColor', 'k', 'LineWidth', 1);

                % plot CI
                CI_Mus = prctile(Mus_data, [2.5, 97.5]);
                CI_Sal = prctile(Sal_data, [2.5, 97.5]);
                CI_Mus_center = mean(CI_Mus);
                CI_Sal_center = mean(CI_Sal);
                CI_Mus_width = CI_Mus(2) - CI_Mus(1);
                CI_Sal_width = CI_Sal(2) - CI_Sal(1);
                errorbar(area_type_idx-0.2, CI_Mus_center, CI_Mus_width/2, 'k', 'LineStyle', 'none', 'LineWidth', 1);
                errorbar(area_type_idx+0.2, CI_Sal_center, CI_Sal_width/2, 'k', 'LineStyle', 'none', 'LineWidth', 1);

                % significance
                p_Mus_left = sum(Mus_data < 0) / length(Mus_data);
                p_Mus_right = sum(Mus_data > 0) / length(Mus_data);
                p_Mus = min(p_Mus_left, p_Mus_right) * 2;
                p_Sal_left = sum(Sal_data < 0) / length(Sal_data);
                p_Sal_right = sum(Sal_data > 0) / length(Sal_data);
                p_Sal = min(p_Sal_left, p_Sal_right) * 2;

                text_x = area_type_idx - 0.2; text_y = CI_Mus(2);
                if p_Mus < 0.001
                    text(text_x, text_y, '***', 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                elseif p_Mus < 0.01
                    text(text_x, text_y, '**', 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                elseif p_Mus < 0.05
                    text(text_x, text_y, '*', 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                end

                text_x = area_type_idx + 0.2; text_y = CI_Sal(2);
                if p_Sal < 0.001
                    text(text_x, text_y, '***', 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                elseif p_Sal < 0.01
                    text(text_x, text_y, '**', 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                elseif p_Sal < 0.05
                    text(text_x, text_y, '*', 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                end

                p_diff_left = sum(Mus_data - Sal_data < 0) / length(Mus_data);
                p_diff_right = sum(Mus_data - Sal_data > 0) / length(Mus_data);
                p_diff = min(p_diff_left, p_diff_right) * 2;

                text_x = area_type_idx; text_y = max(CI_Mus(2), CI_Sal(2))*1.2;
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
            legend({'Muscimol', 'Saline'}, 'Location', 'northwest', 'FontSize', 12, 'Box', 'off');
            ylabel('Value');
            set(gca, 'XTick', 1:n_area_types, 'XTickLabel', area_type_names, 'FontSize', 12);
            xlim([0.5, n_area_types+0.5]);
        end

        % save figure
        file_name = [file_folder, 'Stability_hist_', state, '_', num2str(kernel_idx), '.png'];
        saveas(f, file_name);
        close(f);
    end
end

%% Plot significant count data
f = figure('Position', [100, 100, 1200, 400*n_area_types], 'Visible', 'off');
t = tiledlayout(n_area_types, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
% t = tiledlayout(n_area_types, 3);
title(t, 'Significant Count Data');

for area_type_idx = 1:n_area_types
    for kernel_idx = 1:n_conn_kernel
        nexttile;
        area_type = area_type_names{area_type_idx};
        overlap_index = zeros(3, 2); % (state, session type)
        index_error = zeros(3, 2); % (state, session type)
        for state_idx = 1:n_states
            for session_type_idx = 1:2
                % calculate overlap index: cohen's kappa
                sigcount_mat = significant_count(:, :, area_type_idx, kernel_idx, state_idx, session_type_idx);
                N = sum(sigcount_mat(:));
                margin_x = sum(sigcount_mat, 2);
                margin_y = sum(sigcount_mat, 1); 
                po = (sigcount_mat(1, 1) + sigcount_mat(2, 2) + sigcount_mat(3, 3)) / N; % observed agreement
                pe = sum(margin_x .* margin_y') / (N^2); % expected agreement
                kappa = (po - pe) / (1 - pe);
                overlap_index(state_idx, session_type_idx) = kappa;

                % simplified error calculation
                index_error(state_idx, session_type_idx) = sqrt(po*(1 - po) / (N*(1 - pe)^2));

                % % full error calculation
                % off_diag = 0; % off-diagonal term in error calculation
                % on_diag = 0; % on-diagonal term in error calculation
                % for i = 1:3
                %     for j = 1:3
                %         if i ~= j
                %             off_diag = off_diag + (sigcount_mat(i, j)/N)*((margin_x(i)/N)+(margin_y(j)/N))^2;
                %         end
                %     end
                %     on_diag = on_diag + (sigcount_mat(i, i)/N)*((margin_x(i)/N) + (margin_y(i)/N))^2;
                % end
                % err = sqrt(po*(1-pe)^2 + (1-po)^2*off_diag - 2*(1-po)*(1-pe)*on_diag - (po*pe-2*pe+po)^2)/((1-pe)^2 * sqrt(N));
                % index_error(state_idx, session_type_idx) = err;

                % % full error, version 2
                % A = 0;
                % B = 0;
                % for i = 1:3
                %     for j = 1:3
                %         if i ~= j
                %             B = B + (1-kappa)^2 * (sigcount_mat(i, j) / N) * ((margin_x(i) / N) + (margin_y(j) / N))^2;
                %         end
                %     end
                %     A = A + (sigcount_mat(i, i) / N) * (1-((margin_x(i) / N) + (margin_y(i) / N))*(1-kappa))^2;
                % end
                % C = (kappa - pe*(1-kappa))^2;
                % index_error(state_idx, session_type_idx) = sqrt(A + B - C) / ((1 - pe)^2 * sqrt(N));
            end
        end
        % Bar plot
        bar(overlap_index, 'grouped');
        hold on;
        % Error bars
        ngroup = size(overlap_index, 1);
        nbar = size(overlap_index, 2);
        groupwidth = min(0.8, nbar/(nbar + 1.5));
        for i = 1:nbar
            x = (1:ngroup) - groupwidth/2 + (2*i-1) * groupwidth/(2*nbar);
            errorbar(x, overlap_index(:, i), index_error(:, i), 'k', 'LineStyle', 'none', 'LineWidth', 1);
        end
        hold off;
        title([area_type, ', Kernel ', num2str(kernel_idx)]);
        ylabel('Cohen''s Kappa');
        set(gca, 'XTick', 1:n_states, 'XTickLabel', states, 'FontSize', 8);
        xtickangle(0);
        legend({'Muscimol', 'Saline'}, 'Location', 'northwest', 'FontSize', 12, 'Box', 'off');
        ylim([-0.1, 1]);
    end
end

% save figure
file_folder = [root_path, 'figures/GLM/bootstrap/'];
if ~exist(file_folder, 'dir')
    mkdir(file_folder);
end
file_name = [file_folder, 'Stability_significant_count.png'];
saveas(f, file_name);
close(f);
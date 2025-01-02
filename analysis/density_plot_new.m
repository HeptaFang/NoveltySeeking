session_type = 'Muscimol';
root_path = '../';
kernel = 'DeltaPure';
reg = 'L2=2';
epoch = '2500';
use_filter = true;
filter_threshold = 1;

area_names = {'ACC', 'Thalamus', 'VLPFC'};
% load data
% states = {'Offer1', 'Offer2', 'Decision', 'InfoAnti', 'InfoResp', 'Reward', 'RandomA', 'RandomB'};
% states = {'Offer1', 'Offer2', 'Decision', 'RandomShort', 'RandomLong', 'RandomA', 'RandomB', 'RestOpen', 'RestClose'};
% states = {'RandomA', 'RandomShort', 'RandomLong'};
% states = {'Simulated', 'Simulated_higher', 'RandomLong'};
states = {'Task', 'RestOpen', 'RestClose'};
state_names = {'Task', 'Eyes\_Open', 'Eyes\_Closed'};
aligns = {'AlignFirst', 'AlignLast', 'AlignRandom'};
n_states = length(states);
if strcmp(session_type, 'Muscimol')
    n_session = 10;
else
    n_session = 5;
end

% % load metadata
% if strcmp(prepost, 'Pre')
%     session_stage_full = [session_type, 'Pre', states{1}, '_full'];
% else
%     session_stage_full = [session_type, 'Post', states{1}, '_cortex'];
% end
% 
% % load model
% file_path = [root_path, 'GLM_model/', session_stage_full,...
%         '/GLM_', session_stage_full, '_', '1', '_',...
%         kernel, '_0_', reg, '_', epoch, '.mat'];
% load(file_path,"n_conn_kernel", "kernel_len", "N");
n_conn_kernel = 3;
J_data = cell(3, 3, n_conn_kernel, n_states, 3, n_session, 3); % (area i, area j, kernel, state, align, session, prepost), each cell is a n_area_i x n_area_j matrix
J_data_err = cell(3, 3, n_conn_kernel, n_states, 3, n_session, 3);

all_J = [];
J_state = [];
J_session = [];
J_kernel = [];
J_i = [];
J_j = [];

prepost_all = {'Pre', 'Post', 'Pre'}; % pre cortex, post cortex, pre full
for prepost_idx = 1:3
    prepost = prepost_all{prepost_idx};
    for session_idx = 1:n_session
        fprintf('Loading session %d\n', session_idx);
        for state_idx = 1:n_states
            state = states{state_idx};
            switch prepost_idx
                case 1
                    session_stage = [session_type, 'Pre', state, '_cortex'];
                case 2
                    session_stage = [session_type, 'Post', state, '_cortex'];
                case 3
                    session_stage = [session_type, 'Pre', state, '_full'];
            end

            for align_idx = 1:3
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
                            J_data{i_eff, j_eff, kernel_idx, state_idx, align_idx, session_idx, prepost_idx} = J_area(:);
                            J_data_err{i_eff, j_eff, kernel_idx, state_idx, align_idx, session_idx, prepost_idx} = J_area_err(:);
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
for kernel_idx = 1:n_conn_kernel

    %% data for density plot
    all_data = cell(3, 3, 3, n_states, 3); % (i, j, prepost, state, align)
    all_error = cell(3, 3, 3, n_states, 3); % (i, j, prepost, state, align)

    for i = 1:n_area
        for j = 1:n_area
            % calculate filter: not nan in all states
            filter = cell(n_session, 1);
            for session_idx = 1:n_session
                nan_filter = true(size(J_data{i, j, kernel_idx, 1, 1, session_idx, 3}));
                % significant_filter = false(size(J_data{i, j, kernel_idx, 1, 1, session_idx, 1}));
                for prepost_idx = 1:3
                    if prepost_idx <= 2 && (i == 2 || j == 2) % skip cortex sessions for Thalamus
                        continue;
                    end
                    for state_idx = 1:n_states
                        for align_idx = 1:3
                            data_mat = J_data{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx};
                            error_mat = J_data_err{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx};
                            nan_filter = nan_filter & ~isnan(data_mat) & ~isnan(error_mat);
                            % significant_filter = significant_filter | abs(data_mat) > 2*error_mat; % criterion
                        end
                    end
                end
                filter{session_idx} = nan_filter;
            end

            % store data for density plot
            for prepost_idx = 1:3
                if prepost_idx <= 2 && (i == 2 || j == 2) % skip post sessions for Thalamus
                    continue;
                end
                for state_idx = 1:n_states
                    for align_idx = 1:3
                        data = [];
                        error = [];
                        for session_idx = 1:n_session
                            filter_session = filter{session_idx};
                            data_session = J_data{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx}(filter_session);
                            error_session = J_data_err{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx}(filter_session);
                            data = [data; data_session(:)];
                            error = [error; error_session(:)];
                        end
                        all_data{i, j, prepost_idx, state_idx, align_idx} = data;
                        all_error{i, j, prepost_idx, state_idx, align_idx} = error;
                    end
                end
            end
        end
    end

    %% density plot
    % cmap = brewermap(256, 'Reds');
    nColors = 256;
    cmap = [linspace(1, 1, nColors)', linspace(1, 0, nColors)', linspace(1, 0, nColors)'];

    %% pre vs post
    for state_idx = 1:n_states
        for align_idx = 1:3
            f = figure("Visible", "off", "Position", [0, 0, 800, 800]);
            t = tiledlayout(6, 6, 'TileSpacing', 'compact', 'Padding', 'compact');

            max_log_density = -inf;
            
            for i = 1:n_area
                for j = 1:n_area
                    if i == 2 || j == 2
                        continue;
                    end
                    data_pre = all_data{i, j, 1, state_idx, align_idx};
                    data_post = all_data{i, j, 2, state_idx, align_idx};
                    err_pre = all_error{i, j, 1, state_idx, align_idx};
                    err_post = all_error{i, j, 2, state_idx, align_idx};
                    assert(length(data_pre) == length(data_post), 'Data length mismatch');

                    if use_filter
                        filter = abs(data_pre) > filter_threshold*err_pre | abs(data_post) > filter_threshold*err_post;
                        data_pre = data_pre(filter);
                        data_post = data_post(filter);
                    end

                    % [im, x_edges, y_edges] = histcounts2(data_pre, data_post, plot_edges, plot_edges);
                    % im_log = log1p(im);
                    % imagesc(x_edges, y_edges, im_log');
                    % set(gca, 'YDir', 'normal'); % flip y axis
                    % xlabel('Pre');
                    % ylabel('Post');
                    % colormap(cmap);
                    % title([area_names{i}, ' to ', area_names{j}]);
                    % max_count = max(max_count, max(im(:)));
                    Xdata = data_pre;
                    Ydata = data_post;

                    % Create a 2D histogram
                    [hist_counts, xedges, yedges] = histcounts2(Xdata, Ydata, 50, 'XBinLimits', [-2, 2], 'YBinLimits', [-2, 2]);
                    hist_log = log1p(hist_counts);  % log1p to avoid log(0)
                    max_log_density = max(max_log_density, max(hist_log(:)));  % Track max log density

                    % Create two 1D histogram for x and y axes
                    hist_x = histcounts(Xdata, xedges);
                    hist_y = histcounts(Ydata, yedges);

                    %% Plot the 2D histogram (span: 2x2, upper right)
                    plot_x = i;
                    plot_y = j;
                    if i == 3
                        plot_x = 2;
                    end
                    if j == 3
                        plot_y = 2;
                    end

                    nexttile((plot_x-1)*18 + (plot_y-1)*3 + 2, [2, 2]);

                    imagesc(xedges, yedges, hist_log');
                    set(gca, 'YDir', 'normal'); % flip y axis
                    colormap(cmap);
                    hold on;
                    plot([-2, 2], [-2, 2], 'k--', 'LineWidth', 1);  % Diagonal line
                    
                    axis equal;
                    xlim([-2 2]);
                    ylim([-2 2]);
                    set(gca, 'XTick', [-2 0 2], 'YTick', [-2 0 2]);
                    
                    % % Add gray lines at zero
                    % line([0, 0], [-2, 2], 'Color', 'black', 'LineWidth', 1);
                    % line([-2, 2], [0, 0], 'Color', 'black', 'LineWidth', 1);
                    hold off;
                    
                    % Labeling
                    % xlabel('J\_pre');
                    % ylabel('J\_post');
                    title([area_names{j} ' to ' area_names{i}]);

                    %% Plot the 1D histograms
                    nexttile((plot_x-1)*18 + (plot_y-1)*3 + 12 + 2, [1, 2]);
                    bar(xedges(1:end-1), hist_x, 'FaceColor', 'k');
                    max_count = max(hist_x);
                    hold on;
                    plot([0, 0], [0, max_count*1.05], 'r--', 'LineWidth', 1);
                    hold off;
                    xlim([0, 2]);
                    set(gca, 'YTick', []);
                    set(gca, 'Xtick', []);
                    set(gca, 'Ydir', 'reverse');
                    xlabel('J\_pre');

                    nexttile((plot_x-1)*18 + (plot_y-1)*3 + 1, [2, 1]);
                    barh(yedges(1:end-1), hist_y, 'FaceColor', 'k');
                    max_count = max(hist_y);
                    hold on;
                    plot([0, max_count*1.05], [0, 0], 'r--', 'LineWidth', 1);
                    hold off;
                    ylim([0, 2]);
                    set(gca, 'XTick', []);
                    set(gca, 'Ytick', []);
                    set(gca, 'Xdir', 'reverse');
                    ylabel('J\_post');
                end
            end

            % set all clim to be the same
            for plot_x = 1:2
                for plot_y = 1:2
                    nexttile((plot_x-1)*18 + (plot_y-1)*3 + 2, [2, 2]);
                    clim([0, max_log_density]);
                end
            end

            % Add a global colorbar
            cb = colorbar;
            cb.Layout.Tile = "east";
            ylabel(cb, 'Density');

            log_ticks = log1p([0, 1, 10, 100, 1000, 10000]);
            set(cb, 'Ticks', log_ticks, 'TickLabels', {'0', '1', '10', '100', '1000', '10000'});

            title_str = ['Kernel ', num2str(kernel_idx), ' ', state_names{state_idx}, ' Pre vs Post'];
            filename = [kernel, '_', num2str(kernel_idx), '_', states{state_idx}, '_', aligns{align_idx}];
            sgtitle(title_str);
            fig_folder = [root_path, 'figures/GLM/', session_type, '/density_plot'];
            check_path(fig_folder);
            saveas(f, [fig_folder, '/', filename, '.png']);

            %% histogram of J differences
            f = figure("Visible", "off", "Position", [0, 0, 800, 800]);
            t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

            for i = 1:n_area
                for j = 1:n_area
                    if i == 2 || j == 2
                        continue;
                    end
                    nexttile;

                    data_pre = all_data{i, j, 1, state_idx, align_idx};
                    data_post = all_data{i, j, 2, state_idx, align_idx};
                    err_pre = all_error{i, j, 1, state_idx, align_idx};
                    err_post = all_error{i, j, 2, state_idx, align_idx};
                    assert(length(data_pre) == length(data_post), 'Data length mismatch');

                    if use_filter
                        filter = abs(data_pre) > filter_threshold*err_pre | abs(data_post) > filter_threshold*err_post;
                        data_pre = data_pre(filter);
                        data_post = data_post(filter);
                    end

                    diff_data = data_pre - data_post;
                    fprintf('%d zeros found\n', sum(diff_data==0));
                    diff_data = diff_data(diff_data~=0);

                    histogram(diff_data, 100, 'BinLimits', [-2, 2]);
                    hold on;
                    max_count = max(histcounts(diff_data, 100, 'BinLimits', [-2, 2]));
                    plot([0, 0], [0, max_count*1.05], 'k--', 'LineWidth', 1);
                    mean_diff = mean(diff_data);
                    plot([mean_diff, mean_diff], [0, max_count*1.05], 'r-', 'LineWidth', 0.3);

                    % % ttest
                    % [~, p] = ttest(diff_data);

                    % signrank
                    p = signrank(diff_data);

                    text(0.5, max_count, ['mean: ', num2str(mean_diff)]);
                    text(0.5, max_count*0.9, ['p: ', num2str(p)]);
                    legend({'', 'J_1 = J_2', 'mean'},'Location','northwest');
                    hold off;

                    xlabel('J_{pre} - J_{post}');
                    ylabel('Count');
                    title([area_names{j}, ' to ', area_names{i}]);
                    xlim([-2, 2]);
                end
            end

            sgtitle(title_str);
            saveas(f, [fig_folder, '/', filename, '_diff_hist.png']);
        end
    end

    prepost_strs = {'Pre_cortex', 'Post', 'Pre'};
    %% task vs eyes open vs eyes closed
    for prepost_idx = 1:3
        prepost = prepost_all{prepost_idx};
        prepost_str = prepost_strs{prepost_idx};

        compare_couples = {[1, 2], [1, 3], [2, 3], [3, 2]}; % task vs eyes open, task vs eyes closed, eyes open vs eyes closed
        for compare_idx = 1:3
            compare = compare_couples{compare_idx};
            for align_idx = 1:3
                if prepost_idx <= 2 
                    f = figure("Visible", "off", "Position", [0, 0, 800, 800]);
                    tiledlayout(6, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
                    row_length = 6;
                    plot_size = 2;
                else
                    f = figure("Visible", "off", "Position", [0, 0, 1200, 1200]);
                    tiledlayout(9, 9, 'TileSpacing', 'compact', 'Padding', 'compact');
                    row_length = 9;
                    plot_size = 3;
                end

                max_log_density = -inf;
                
                for i = 1:n_area
                    for j = 1:n_area
                        if (i == 2 || j == 2) && prepost_idx <= 2 % skip cortex sessions for Thalamus
                            continue;
                        end
                        Xdata = all_data{i, j, prepost_idx, compare(1), align_idx};
                        Ydata = all_data{i, j, prepost_idx, compare(2), align_idx};
                        Xerr = all_error{i, j, prepost_idx, compare(1), align_idx};
                        Yerr = all_error{i, j, prepost_idx, compare(2), align_idx};
                        assert(length(Xdata) == length(Ydata), 'Data length mismatch');

                        if use_filter
                            filter = abs(Xdata) > filter_threshold*Xerr | abs(Ydata) > filter_threshold*Yerr;
                            Xdata = Xdata(filter);
                            Ydata = Ydata(filter);
                        end

                        % Create a 2D histogram
                        [hist_counts, xedges, yedges] = histcounts2(Xdata, Ydata, 50, 'XBinLimits', [-2, 2], 'YBinLimits', [-2, 2]);
                        hist_log = log1p(hist_counts);  % log1p to avoid log(0)
                        max_log_density = max(max_log_density, max(hist_log(:)));  % Track max log density

                        % Create two 1D histogram for x and y axes
                        hist_x = histcounts(Xdata, xedges);
                        hist_y = histcounts(Ydata, yedges);

                        %% Plot the 2D histogram (span: 2x2, upper right)
                        plot_x = i;
                        plot_y = j;
                        if i == 3 && prepost_idx <= 2
                            plot_x = 2;
                        end
                        if j == 3 && prepost_idx <= 2
                            plot_y = 2;
                        end

                        nexttile((plot_x-1)*3*row_length + (plot_y-1)*3 + 2, [2, 2]);

                        imagesc(xedges, yedges, hist_log');
                        set(gca, 'YDir', 'normal'); % flip y axis
                        colormap(cmap);
                        hold on;
                        plot([-2, 2], [-2, 2], 'k--', 'LineWidth', 1);  % Diagonal line
                        
                        axis equal;
                        xlim([-2 2]);
                        ylim([-2 2]);
                        set(gca, 'XTick', [-2 0 2], 'YTick', [-2 0 2]);
                        
                        % Add gray lines at zero
                        % line([0, 0], [-2, 2], 'Color', 'black', 'LineWidth', 1);
                        % line([-2, 2], [0, 0], 'Color', 'black', 'LineWidth', 1);
                        hold off;
                        
                        % Labeling
                        % xlabel('J\_pre');
                        % ylabel('J\_post');
                        title([area_names{j} ' to ' area_names{i}]);

                        %% Plot the 1D histograms
                        nexttile((plot_x-1)*3*row_length + (plot_y-1)*3 + 2*row_length + 2, [1, 2]);
                        bar(xedges(1:end-1), hist_x, 'FaceColor', 'k');
                        max_count = max(hist_x);
                        hold on;
                        plot([0, 0], [0, max_count*1.05], 'r--', 'LineWidth', 1);
                        hold off;
                        xlim([0, 2]);
                        set(gca, 'YTick', []);
                        set(gca, 'Xtick', []);
                        set(gca, 'Ydir', 'reverse');
                        xlabel(states{compare(1)});

                        nexttile((plot_x-1)*3*row_length + (plot_y-1)*3 + 1, [2, 1]);
                        barh(yedges(1:end-1), hist_y, 'FaceColor', 'k');
                        max_count = max(hist_y);
                        hold on;
                        plot([0, max_count*1.05], [0, 0], 'r--', 'LineWidth', 1);
                        hold off;
                        ylim([-2, 2]);
                        set(gca, 'XTick', []);
                        set(gca, 'Ytick', []);
                        set(gca, 'Xdir', 'reverse');
                        ylabel(states{compare(2)});
                    end
                end

                % set all clim to be the same
                for plot_x = 1:plot_size
                    for plot_y = 1:plot_size
                        nexttile((plot_x-1)*3*row_length + (plot_y-1)*3 + 2, [2, 2]);
                        clim([0, max_log_density]);
                    end
                end

                % Add a global colorbar
                cb = colorbar;
                cb.Layout.Tile = "east";
                ylabel(cb, 'Density');

                log_ticks = log1p([0, 1, 10, 100, 1000, 10000]);
                set(cb, 'Ticks', log_ticks, 'TickLabels', {'0', '1', '10', '100', '1000', '10000'});

                title_str = ['Kernel ', num2str(kernel_idx), ' ', prepost_str, ' ', state_names{compare(1)}, ' vs ', state_names{compare(2)}];
                filename = [kernel, '_', num2str(kernel_idx), '_', prepost_str, '_', states{compare(1)}, '_', states{compare(2)}, '_', aligns{align_idx}];
                sgtitle(title_str);
                fig_folder = [root_path, 'figures/GLM/', session_type, '/density_plot'];
                check_path(fig_folder);
                saveas(f, [fig_folder, '/', filename, '.png']);

                %% histogram of J differences
                if prepost_idx <= 2 
                    f = figure("Visible", "off", "Position", [0, 0, 800, 800]);
                    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
                else
                    f = figure("Visible", "off", "Position", [0, 0, 1200, 1200]);
                    tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
                end

                for i = 1:n_area
                    for j = 1:n_area
                        if (i == 2 || j == 2) && prepost_idx <= 2 % skip cortex sessions for Thalamus
                            continue;
                        end
                        nexttile;

                        Xdata = all_data{i, j, prepost_idx, compare(1), align_idx};
                        Ydata = all_data{i, j, prepost_idx, compare(2), align_idx};
                        Xerr = all_error{i, j, prepost_idx, compare(1), align_idx};
                        Yerr = all_error{i, j, prepost_idx, compare(2), align_idx};
                        assert(length(Xdata) == length(Ydata), 'Data length mismatch');

                        if use_filter
                            filter = abs(Xdata) > filter_threshold*Xerr | abs(Ydata) > filter_threshold*Yerr;
                            Xdata = Xdata(filter);
                            Ydata = Ydata(filter);
                        end

                        diff_data = Xdata - Ydata;
                        fprintf('%d zeros found\n', sum(diff_data==0));
                        diff_data = diff_data(diff_data~=0);

                        histogram(diff_data, 100, 'BinLimits', [-2, 2]);
                        hold on;
                        max_count = max(histcounts(diff_data, 100, 'BinLimits', [-2, 2]));
                        plot([0, 0], [0, max_count*1.05], 'k--', 'LineWidth', 1);
                        mean_diff = mean(diff_data);
                        plot([mean_diff, mean_diff], [0, max_count*1.05], 'r-', 'LineWidth', 0.3);

                        % ttest
                        % [~, p] = ttest(diff_data);

                        % signrank
                        p = signrank(diff_data);

                        text(0.5, max_count, ['mean: ', num2str(mean_diff)]);
                        text(0.5, max_count*0.9, ['p: ', num2str(p)]);

                        hold off;

                        xlabel([states{compare(1)}, ' - ', states{compare(2)}]);
                        ylabel('Count');
                        title([area_names{j}, ' to ', area_names{i}]);
                        xlim([-2, 2]);
                        legend({'', 'J_1 = J_2', 'mean'},'Location','northwest');
                    end
                end

                sgtitle(title_str);
                saveas(f, [fig_folder, '/', filename, '_diff_hist.png']);
            end
        end
    end

end
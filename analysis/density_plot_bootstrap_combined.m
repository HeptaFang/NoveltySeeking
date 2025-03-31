% Combine areawise data to across-area/within-area data for density plot

% Todo: Mode=Mean, fix across-area data.
% Todo: Plotting logic: Musimol within vs Saline within, Musimol across vs Saline across.

root_path = '../';
kernel = 'DeltaPure';
reg = 'L2=2';
epoch = '2500';
use_filter =false;
filter_threshold = 1;
plot_mode = 'Scatter'; % 'Scatter' or 'Density'
% plot_mode = 'Density'; % 'Scatter' or 'Density'
data_mode = 'Individual'; % 'Individual' or 'Mean'
% data_mode = 'Mean'; % 'Individual' or 'Mean'

% plot_range = [-3, 3];
plot_range_within = [-3, 3];
plot_range_across = [-1.5, 1.5];

addpath("UtilFuncs");
addpath("UtilFuncs/HELPER_GENERAL");

area_names = {'ACC', 'Thalamus', 'VLPFC'};
% load data
% states = {'Offer1', 'Offer2', 'Decision', 'InfoAnti', 'InfoResp', 'Reward', 'RandomA', 'RandomB'};
% states = {'Offer1', 'Offer2', 'Decision', 'RandomShort', 'RandomLong', 'RandomA', 'RandomB', 'RestOpen', 'RestClose'};
% states = {'RandomA', 'RandomShort', 'RandomLong'};
% states = {'Simulated', 'Simulated_higher', 'RandomLong'};
states = {'Task', 'RestOpen', 'RestClose'};
state_names = {'Task', 'Eyes\_Open', 'Eyes\_Closed'};
session_types = {'Muscimol', 'Saline'};
% aligns = {'AlignFirst', 'AlignLast', 'AlignRandom'};
aligns = {'AlignLast'};
n_states = length(states);
n_aligns = length(aligns);

n_conn_kernel = 3;
J_data = cell(3, 3, n_conn_kernel, n_states, n_aligns, 10, 3, 2); % (area i, area j, kernel, state, align, session, prepost, session_type), each cell is a n_area_i x n_area_j matrix
J_data_err = cell(3, 3, n_conn_kernel, n_states, n_aligns, 10, 3, 2);

for session_type_idx = 1:2
    session_type = session_types{session_type_idx};
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

    % all_J = [];
    % J_state = [];
    % J_session = [];
    % J_kernel = [];
    % J_i = [];
    % J_j = [];

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
                                J_data{i_eff, j_eff, kernel_idx, state_idx, align_idx, session_idx, prepost_idx, session_type_idx} = J_area;
                                J_data_err{i_eff, j_eff, kernel_idx, state_idx, align_idx, session_idx, prepost_idx, session_type_idx} = J_area_err;
                                % all_J = [all_J; J_area(:)];
                                % J_state = [J_state, repelem(state_idx, numel(J_area))];
                                % J_session = [J_session, repelem(session_idx, numel(J_area))];
                                % J_kernel = [J_kernel, repelem(kernel_idx, numel(J_area))];
                                % J_i = [J_i, repelem(i_eff, numel(J_area))];
                                % J_j = [J_j, repelem(j_eff, numel(J_area))];
                            end
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
    all_data = cell(2, 3, n_states, 3, 2); % (within/across, prepost, state, align, session_type)
    all_error = cell(2, 3, n_states, 3, 2); % (within/across, prepost, state, align, session_type)
    for session_type_idx = 1:2
        session_type = session_types{session_type_idx};
        if strcmp(session_type, 'Muscimol')
            n_session = 10;
        else
            n_session = 5;
        end
        for i = 1:n_area
            for j = 1:n_area
                % calculate filter: not nan in all states
                filter = cell(n_session, 1);
                for session_idx = 1:n_session
                    nan_filter = true(size(J_data{i, j, kernel_idx, 1, 1, session_idx, 3, session_type_idx})); % all nan
                    % significant_filter = false(size(J_data{i, j, kernel_idx, 1, 1, session_idx, 1}));
                    for prepost_idx = 1:3
                        if prepost_idx <= 2 && (i == 2 || j == 2) % skip cortex sessions for Thalamus
                            continue;
                        end
                        for state_idx = 1:n_states
                            for align_idx = 1:n_aligns
                                data_mat = J_data{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx, session_type_idx};
                                error_mat = J_data_err{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx, session_type_idx};
                                nan_filter = nan_filter & ~isnan(data_mat) & ~isnan(error_mat);
                                % significant_filter = significant_filter | abs(data_mat) > 2*error_mat; % criterion
                            end
                        end
                    end
                    filter{session_idx} = nan_filter;
                end

                % store data for scatter/density plot
                for prepost_idx = 1:3
                    if prepost_idx <= 2 && (i == 2 || j == 2) % skip post sessions for Thalamus
                        continue;
                    end
                    for state_idx = 1:n_states
                        for align_idx = 1:n_aligns
                            data = [];
                            error = [];
                            for session_idx = 1:n_session
                                filter_session = filter{session_idx};
                                switch data_mode
                                    case 'Individual'
                                        % use individual Js
                                        data_session = J_data{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx, session_type_idx}(filter_session);
                                        error_session = J_data_err{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx, session_type_idx}(filter_session);
                                    case 'Mean'
                                        % use Js meant by cell
                                        J_mat = J_data{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx, session_type_idx};
                                        J_err_mat = J_data_err{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx, session_type_idx};
                                        J_in = sum(J_mat, 2, 'omitnan')/size(J_mat, 1);
                                        J_out = sum(J_mat, 1, 'omitnan')/size(J_mat, 2);

                                        % This error calculation is not correct. To be fixed.
                                        J_err_in = mean(J_err_mat, 1, 'omitnan');
                                        J_err_out = mean(J_err_mat, 2, 'omitnan');
                                        
                                        if i==j
                                            data_session = (J_in(:) + J_out(:)) / 2;
                                        else
                                            data_session = J_in(:);
                                        end

                                        % error_session = sqrt(J_err_in(:).^2 + J_err_out(:).^2) / 2;
                                        error_session = J_err_in(:);
                                end
                                data = [data; data_session(:)];
                                error = [error; error_session(:)];
                            end
                            if i==j
                                all_data{1, prepost_idx, state_idx, align_idx, session_type_idx} = [all_data{1, prepost_idx, state_idx, align_idx, session_type_idx};data];
                                all_error{1, prepost_idx, state_idx, align_idx, session_type_idx} = [all_error{1, prepost_idx, state_idx, align_idx, session_type_idx};error];
                            else
                                all_data{2, prepost_idx, state_idx, align_idx, session_type_idx} = [all_data{2, prepost_idx, state_idx, align_idx, session_type_idx};data];
                                all_error{2, prepost_idx, state_idx, align_idx, session_type_idx} = [all_error{2, prepost_idx, state_idx, align_idx, session_type_idx};error];
                            end
                        end
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
        for align_idx = 1:n_aligns
            f = figure("Visible", "off", 'PaperPosition', [0 0 15 15]);
            t = tiledlayout(6, 6, 'TileSpacing', 'compact', 'Padding', 'compact');

            max_log_density = -inf;
            
            for i = 1:n_area
                for j = 1:n_area
                    if i == 2 || j == 2
                        continue;
                    end

                    bootstrap_data = cell(2, 2);

                    for session_type_idx = 1:2
                        data_pre = all_data{i, j, 1, state_idx, align_idx, session_type_idx};
                        data_post = all_data{i, j, 2, state_idx, align_idx, session_type_idx};
                        err_pre = all_error{i, j, 1, state_idx, align_idx, session_type_idx};
                        err_post = all_error{i, j, 2, state_idx, align_idx, session_type_idx};
                        assert(length(data_pre) == length(data_post), 'Data length mismatch');

                        if use_filter
                            filter = abs(data_pre) > filter_threshold*err_pre | abs(data_post) > filter_threshold*err_post;
                            data_pre = data_pre(filter);
                            data_post = data_post(filter);
                        end
                        Xdata = data_pre;
                        Ydata = data_post;

                        bootstrap_data{session_type_idx, 1} = Xdata;
                        bootstrap_data{session_type_idx, 2} = Ydata;
                        %% Plot the 2D histogram (span: 2x2, upper right)
                        plot_x = i;
                        plot_y = j;
                        if i == 3
                            plot_x = 2;
                        end
                        if j == 3
                            plot_y = 2;
                        end
                        nexttile((plot_x-1)*36 + (plot_y-1)*6 + 2 + 3*(session_type_idx-1), [2, 2]);

                        if i==j
                            plot_range = plot_range_within;
                        else
                            plot_range = plot_range_across;
                        end

                        switch plot_mode
                            case 'Density'
                                % Create a 2D histogram
                                [hist_counts, Xedges, Yedges] = histcounts2(Xdata, Ydata, 100, 'XBinLimits', plot_range, 'YBinLimits', plot_range);
                                hist_log = log1p(hist_counts);  % log1p to avoid log(0)
                                max_log_density = max(max_log_density, max(hist_log(:)));  % Track max log density

                                imagesc(xedges, yedges, hist_log');
                                colormap(cmap);
                                set(gca, 'YDir', 'normal'); % flip y axis

                            case 'Scatter'
                                scatter(Xdata, Ydata, 15, '+');
                                % scatter(Xdata, Ydata, 5, '+');
                                xedges = -2:0.04:2;
                                yedges = -2:0.04:2;

                        end

                        hold on;
                        plot(plot_range, plot_range, 'k--', 'LineWidth', 1);  % Diagonal line
                        
                        axis equal;
                        xlim(plot_range);
                        ylim(plot_range);
                        set(gca, 'XTick', [plot_range(1) 0 plot_range(2)], 'YTick', [plot_range(1) 0 plot_range(2)]);
                        
                        % Add gray lines at zero
                        line([0, 0], plot_range, 'Color', 'black', 'LineWidth', 1);
                        line(plot_range, [0, 0], 'Color', 'black', 'LineWidth', 1);
                        hold off;
                        
                        
                        % Labeling
                        % xlabel('J\_pre');
                        % ylabel('J\_post');
                        title([area_names{j} ' to ' area_names{i}, ', ', session_types{session_type_idx}]);

                        %% Plot the 1D histograms
                        % Create two 1D histogram for x and y axes
                        hist_x = histcounts(Xdata, xedges);
                        hist_y = histcounts(Ydata, yedges);

                        % plot
                        nexttile((plot_x-1)*36 + (plot_y-1)*6 + 24 + 2 + 3*(session_type_idx-1), [1, 2]);
                        bar(xedges(1:end-1), hist_x, 'FaceColor', 'k');
                        max_count = max(hist_x);
                        hold on;
                        plot([0, 0], [0, max_count*1.05], 'r--', 'LineWidth', 1);
                        hold off;
                        xlim(plot_range);
                        set(gca, 'YTick', []);
                        set(gca, 'Xtick', []);
                        set(gca, 'Ydir', 'reverse');
                        xlabel('J\_pre');

                        nexttile((plot_x-1)*36 + (plot_y-1)*6 + 1 + 3*(session_type_idx-1), [2, 1]);
                        barh(yedges(1:end-1), hist_y, 'FaceColor', 'k');
                        max_count = max(hist_y);
                        hold on;
                        plot([0, max_count*1.05], [0, 0], 'r--', 'LineWidth', 1);
                        hold off;
                        ylim(plot_range);
                        set(gca, 'XTick', []);
                        set(gca, 'Ytick', []);
                        set(gca, 'Xdir', 'reverse');
                        ylabel('J\_post');
                    end

                    % Bootstrap
                    [p, ci, ci1, ci2] = CorrCompareInd(bootstrap_data{1, 1}, bootstrap_data{1, 2}, bootstrap_data{2, 1}, bootstrap_data{2, 2}, 1000);
                    % plot on figure
                    nexttile((plot_x-1)*36 + (plot_y-1)*6 + 24+4);
                    text(0.5, 0.7, ['p = ', num2str(p)], 'FontSize', 10, 'HorizontalAlignment', 'center');
                    text(0.5, 0.6, ['CI, diff: [', num2str(ci(1)), ', ', num2str(ci(2)), ']'], 'FontSize', 10, 'HorizontalAlignment', 'center');
                    text(0.5, 0.5, ['CI, Muscimol: [', num2str(ci1(1)), ', ', num2str(ci1(2)), ']'], 'FontSize', 10, 'HorizontalAlignment', 'center');
                    text(0.5, 0.4, ['CI, Saline: [', num2str(ci2(1)), ', ', num2str(ci2(2)), ']'], 'FontSize', 10, 'HorizontalAlignment', 'center');
                    axis off;
                end
            end

            % set all clim to be the same
            % for plot_x = 1:2
            %     for plot_y = 1:2
            %         nexttile((plot_x-1)*18 + (plot_y-1)*3 + 2, [2, 2]);
            %         clim([0, max_log_density]);
            %     end
            % end

            % Add a global colorbar
            % cb = colorbar;
            % cb.Layout.Tile = "east";
            % ylabel(cb, 'Density');

            % log_ticks = log1p([0, 1, 10, 100, 1000, 10000]);
            % set(cb, 'Ticks', log_ticks, 'TickLabels', {'0', '1', '10', '100', '1000', '10000'});

            title_str = ['Kernel ', num2str(kernel_idx), ' ', state_names{state_idx}, ' Pre vs Post'];
            filename = [kernel, '_', num2str(kernel_idx), '_', states{state_idx}, '_', aligns{align_idx}];
            sgtitle(title_str);
            fig_folder = [root_path, 'figures/GLM/bootstrap/density_plot'];
            check_path(fig_folder);
            % saveas(f, [fig_folder, '/', filename, '.png']);
            print(f, [fig_folder, '/', filename, '.png'], '-dpng', '-r300');
        end
    end
end
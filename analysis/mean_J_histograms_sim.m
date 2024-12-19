root_path = '../';
kernel = 'DeltaPure';
epoch = '2500';
states = {'Simulated_lower', 'Simulated', 'Simulated_higher'};
n_states = length(states);
n_session = 1;
session_type = 'Simulated';
addpath("figtool/");

n_conn_kernel = 3;
regs = {'L1=2', 'L2=2', 'NoReg'};
J_data = cell(n_conn_kernel, n_states, length(regs)+1); % (kernel, state, reg/noreg/actual), each cell is a n_area_i x n_area_j matrix
J_data_err = cell(n_conn_kernel, n_states, length(regs)+1);
filters = cell(n_conn_kernel);
firing_rate_data = cell(n_states);

all_J = [];
J_state = [];
J_session = [];
J_kernel = [];
J_i = [];
J_j = [];

%% load data to J_data and J_data_err
for session_idx = 1:n_session
    fprintf('Loading session %d\n', session_idx);
    for state_idx = 1:n_states
        state = states{state_idx};
        session_stage_full = state;

        %% load reg and noreg data
        for reg_idx = 1:length(regs)
            reg = regs{reg_idx};

            file_path = [root_path, 'GLM_model/', session_stage_full,...
                '/GLM_', session_stage_full, '_', num2str(session_idx), '_',...
                kernel, '_0_', reg, '_', epoch, '.mat'];

            fprintf('Loading %s\n', session_stage_full);

            load(file_path, "model_par", "n_PS_kernel", "kernel_len", "N", "model_err");

            for kernel_idx = 1:n_conn_kernel
                J_mat = model_par(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));

                % temp use: fix this
                if isa(model_err, 'struct')
                    % model_err = model_err.minuslogL;
                    model_err = model_err.total;
                end

                J_err = model_err(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));

                J_data{kernel_idx, state_idx, reg_idx} = J_mat;
                J_data_err{kernel_idx, state_idx, reg_idx} = J_err;
            end
        end

        %% load actual data
        file_path = [root_path, 'GLM_data/', state, '/raster_', state, '_1_0.mat'];
        J_actual = load(file_path).J;
        for kernel_idx = 1:n_conn_kernel
            J_data{kernel_idx, state_idx, length(regs)+1} = J_actual(:, :, kernel_idx);
        end

        %% load firing rate data
        file_path = [root_path, 'GLM_data/', state, '/GLMdata_', state, '_1_', kernel, '_0.mat'];
        r = load(file_path).raster;
        firing_rate_data{state_idx} = mean(r, 2)*1000;
    end
end

%% plot J histograms
ylim_all_kernel = {[-0.15, 1], [-0.01, 0.15], [-0.01, 0.12]};
for kernel_idx = 1:n_conn_kernel
    % Each kernel is one figure
    f = figure("Visible", "off","Position",[0, 0, 1600, 600]);

    %% J panel
    subplot('Position', [0.05, 0.1, 0.6, 0.8]);
    title('Average J with error bars');
    
    % barplot
    means = zeros(n_states, length(regs)+1); % (state, reg/noreg/actual)
    errors = zeros(n_states, length(regs)+1);

    % calculate filter: not nan in all states and significant in at least one state
    nan_filter = true(size(J_data{kernel_idx, 1, 1}));
    significant_filter = false(size(J_data{kernel_idx, 1, 1}));
    for reg_idx = 1:2
        for state_idx = 1:n_states
            data_mat = J_data{kernel_idx, state_idx, reg_idx};
            error_mat = J_data_err{kernel_idx, state_idx, reg_idx};
            nan_filter = nan_filter & ~isnan(data_mat) & ~isnan(error_mat);
            significant_filter = significant_filter | abs(data_mat) > 2*error_mat; % criterion
        end
    end
    filter = nan_filter & significant_filter;
    filters{kernel_idx} = filter;

    % calculate mean and error 
    for reg_idx = 1:(length(regs)+1)
        for state_idx = 1:n_states
            data = J_data{kernel_idx, state_idx, reg_idx}(filter);

            % nan_filter = ~isnan(data) & ~isnan(error);
            % significant_filter = abs(data) > 2*error;
            % % filter = nan_filter & significant_filter;
            % filter = nan_filter;

            J_mean = mean(data, "all");
            
            if reg_idx < length(regs)+1
                error = J_data_err{kernel_idx, state_idx, reg_idx}(filter);
                % Variance of the data points (sample standard deviation)
                data_variance = var(data, 0, "all");
                sem_data = sqrt(data_variance / numel(data)); % Standard error of the mean

                % Propagated error from individual measurement errors
                propagated_error = sqrt(sum(error.^2) / numel(data)^2);

                % Combine the two errors
                combined_error = sqrt(sem_data^2 + propagated_error^2);

                means(state_idx, reg_idx) = J_mean;
                errors(state_idx, reg_idx) = combined_error;
            else
                % actual data
                data_variance = var(data, 0, "all");
                sem_data = sqrt(data_variance / numel(data)); % Standard error of the mean

                means(state_idx, reg_idx) = J_mean;
                errors(state_idx, reg_idx) = sem_data;
            end
        end
    end

    bar_width = 0.25;
    for state_idx = 1:n_states
        % distribute evenly
        shift = (state_idx - (n_states+1)/2) * bar_width;
        hold on;
        bar((1:(length(regs)+1))+shift, means(state_idx, :), bar_width, '');
        hold off;
    end
    for state_idx = 1:n_states
        % statistical test with the actual data
        % distribute evenly
        shift = (state_idx - (n_states+1)/2) * bar_width;
        hold on;
        errorbar((1:(length(regs)+1))+shift, means(state_idx, :), errors(state_idx, :), 'LineStyle', 'none', "Color", [0 0 0],"CapSize",2);

        for reg_idx = 1:length(regs)
            data = J_data{kernel_idx, state_idx, reg_idx}(filter);
            actual_data = J_data{kernel_idx, state_idx, length(regs)+1}(filter);
            [~, p] = ttest(data, actual_data);
            % p = signrank(data, actual_data);
            if p < 0.001
                text(reg_idx+shift, means(state_idx, reg_idx)+errors(state_idx, reg_idx)+0.01, '***', 'HorizontalAlignment', 'center');
            elseif p < 0.01
                text(reg_idx+shift, means(state_idx, reg_idx)+errors(state_idx, reg_idx)+0.01, '**', 'HorizontalAlignment', 'center');
            elseif p < 0.05
                text(reg_idx+shift, means(state_idx, reg_idx)+errors(state_idx, reg_idx)+0.01, '*', 'HorizontalAlignment', 'center');
            end

            if state_idx < n_states
                % test adjacent states
                data = J_data{kernel_idx, state_idx, reg_idx}(filter);
                data_next = J_data{kernel_idx, state_idx+1, reg_idx}(filter);
                [~, p] = ttest(data, data_next);
                % p = signrank(data, data_next);
                line_height = max([means(state_idx, reg_idx)+errors(state_idx, reg_idx), means(state_idx+1, reg_idx)+errors(state_idx+1, reg_idx)]) + 0.05;
                if p < 0.001
                    % horizontal line and stars at the top
                    plot([reg_idx+shift, reg_idx+shift+bar_width], [line_height, line_height], 'k-');
                    text(reg_idx+shift+bar_width/2, line_height+0.02, '***', 'HorizontalAlignment', 'center');
                elseif p < 0.01
                    plot([reg_idx+shift, reg_idx+shift+bar_width], [line_height, line_height], 'k-');
                    text(reg_idx+shift+bar_width/2, line_height+0.02, '**', 'HorizontalAlignment', 'center');
                elseif p < 0.05
                    plot([reg_idx+shift, reg_idx+shift+bar_width], [line_height, line_height], 'k-');
                    text(reg_idx+shift+bar_width/2, line_height+0.02, '*', 'HorizontalAlignment', 'center');
                end
            end

        end

        hold off;
    end
    title(['Average J for kernel ', num2str(kernel_idx)]);
    xticks(1:(length(regs)+1));
    xticklabels(["L1 regularization", "L2 regularization", "No regularization", "Actual J"]);
    xlabel('Model regularization type')
    ylabel('Average J');

    legend_strings = {'Low firing rate dataset', 'Mid firing rate dataset', 'High firing rate dataset'};
    legend(legend_strings); 
    % ylim(ylim_all_kernel{kernel_idx});
    ylim([0, 0.65]);

    % Create the right panel
    subplot('Position', [0.05 + 0.6 + 0.1, 0.1, 0.2, 0.8]);

    bar_width = 0.25;
    for state_idx = 1:n_states
        firing_rate = firing_rate_data{state_idx};
        mean_firing_rate = mean(firing_rate);
        err_firing_rate = std(firing_rate) / sqrt(numel(firing_rate));

        % distribute evenly
        shift = (state_idx - (n_states+1)/2) * bar_width;
        hold on;
        bar(1+shift, mean_firing_rate, bar_width, '');
        hold off;
    end
    for state_idx = 1:n_states
        firing_rate = firing_rate_data{state_idx};
        mean_firing_rate = mean(firing_rate);
        err_firing_rate = std(firing_rate) / sqrt(numel(firing_rate));

        % distribute evenly
        shift = (state_idx - (n_states+1)/2) * bar_width;
        hold on;
        errorbar(1+shift, mean_firing_rate, err_firing_rate, 'LineStyle', 'none', "Color", [0 0 0],"CapSize",2);
        hold off;
    end
    title(['Firing rate for kernel ', num2str(kernel_idx)]);
    xticks([]);
    ylabel('Firing rate (Hz)');
    title('Firing rate');
    legend(legend_strings);

    fig_folder = [root_path, 'figures/GLM/', session_type];
    check_path(fig_folder);
    saveas(f, [fig_folder, '/J_histograms_sim_', num2str(kernel_idx), '.png']);
end

%% plot J matrices
cmap = brewermap(256, '*RdBu');
% cmap = brewermap(256, 'RdBu');
% cmap = [0 0 0; cmap];
for kernel_idx = 1:n_conn_kernel
    % 3X6 subplot. rows: low firing rate, mid firing rate, high firing rate
    % columns: reg, noreg, actual, filtered reg, filtered noreg, filtered actual
    f = figure("Visible", "off","Position",[0, 0, 1600, 800]);
    for state_idx = 1:n_states
        for reg_idx = 1:(length(regs)+1)
            subplot(3, 2*(length(regs)+1), (state_idx-1)*2*(length(regs)+1) + reg_idx);
            J = J_data{kernel_idx, state_idx, reg_idx};
            imagesc(J, 'AlphaData', ~isnan(J));
            axis square;
            regs_strings = [regs, 'Actual'];
            if state_idx == 1
                title(['J matrix, ', regs_strings{reg_idx}]);
            end
            if reg_idx == 1
                ylabel([legend_strings{state_idx}]);
            end
            clim([-3, 3]);
            colormap(cmap);
            xticks([]);
            yticks([]);
        end
        for reg_idx = 1:(length(regs)+1)
            subplot(3, 2*(length(regs)+1), (state_idx-1)*2*(length(regs)+1) + reg_idx + length(regs)+1);
            J = J_data{kernel_idx, state_idx, reg_idx};
            J(~filters{kernel_idx}) = NaN;
            % set(gca, 'Color', 'k');
            % set(f, 'Color', 'k');
            imagesc(J, 'AlphaData', ~isnan(J));
            % imagesc(J);
            axis square;
            regs_strings = [regs, 'Actual'];
            if state_idx == 1
                title(['Filtered J, ', regs_strings{reg_idx}]);
            end
            clim([-3, 3]);
            colormap(cmap);
            xticks([]);
            yticks([]);
        end
    end

    sgtitle(['J matrices for kernel ', num2str(kernel_idx)]);
    fig_folder = [root_path, 'figures/GLM/', session_type];
    check_path(fig_folder);
    saveas(f, [fig_folder, '/J_matrices_sim_', num2str(kernel_idx), '.png']);
    % export_fig(f, [fig_folder, '/J_matrices_sim_', num2str(kernel_idx), '.png']);
end
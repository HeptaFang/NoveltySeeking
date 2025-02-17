% Use last aligned time bins
% Use 3x3 grid to show pre-post significant count change
% Plot by session

session_type = 'Muscimol';
root_path = '../';
kernel = 'DeltaPure';
reg = 'L2=2';
epoch = '2500';
area_names = {'ACC', 'Thalamus', 'VLPFC'};
filter_threshold = 1;
% load data
% states = {'Offer1', 'Offer2', 'Decision', 'InfoAnti', 'InfoResp', 'Reward', 'RandomA', 'RandomB'};
% states = {'Offer1', 'Offer2', 'Decision', 'RandomShort', 'RandomLong', 'RandomA', 'RandomB', 'RestOpen', 'RestClose'};
% states = {'RandomA', 'RandomShort', 'RandomLong'};
% states = {'Simulated', 'Simulated_higher', 'RandomLong'};
states = {'Task', 'RestOpen', 'RestClose'};
aligns = {'AlignLast'};
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
J_data = cell(3, 3, n_conn_kernel, n_states, 1, n_session, 2); % (area i, area j, kernel, state, align, session, prepost), each cell is a n_area_i x n_area_j matrix
J_data_err = cell(3, 3, n_conn_kernel, n_states, 1, n_session, 2);

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
ylim_all_kernel = {[-0.05, 0.9], [-0.15, 1], [-0.01, 0.12]};
for kernel_idx = 1:n_conn_kernel
    % Each kernel is one figure
    % f = figure("Visible", "off","Position",[0, 0, 1600, 900]);
    % % ij are subplots
    % tiles = tiledlayout(n_area, n_area);

    % for i = 1:n_area
    %     for j = 1:n_area
    %         fprintf('Kernel%d, i=%d, j=%d\n', kernel_idx, i, j);
    %         nexttile;
    %         % % violin plot
    %         % filter = J_i == i & J_j == j & J_kernel == kernel_idx;
    %         % data = all_J(filter);
    %         % data_state = J_state(filter);
    %         % data_session = J_session(filter);
    %         % data_by_session = cell(1, n_session);
    %         % for session_idx = 1:n_session
    %         %     data_by_session{session_idx} = data(data_session == session_idx);
    %         %     if size(data_by_session{session_idx}, 1) == 0
    %         %         data_by_session{session_idx} = 0;
    %         %     end
    %         % end
    %         % violinplot(data_session, data, GroupByColor=data_state);
    %         % % daviolinplot(data_by_session,'groups',data_state,'outsymbol','k+',...
    %         % %     'xtlabels', 1:n_session,'scatter',2,'jitter',1,...
    %         % %     'box',1,'boxcolors','same', 'color', c,'scattercolors','same',...
    %         % %     'boxspacing',1.1,'legend',states);
    %         % % legend(states);

    %         % barplot
    %         means = zeros(2, n_states, 3);
    %         errors = zeros(2, n_states, 3);

    %         % calculate filter: not nan in all states and significant in at least one state

    %         filter = cell(n_session, 1);
    %         for session_idx = 1:n_session
    %             nan_filter = true(size(J_data{i, j, kernel_idx, 1, 1, session_idx, 1}));
    %             significant_filter = false(size(J_data{i, j, kernel_idx, 1, 1, session_idx, 1}));
    %             for prepost_idx = 1:2
    %                 if prepost_idx == 2 && (i == 2 || j == 2) % skip post sessions for Thalamus
    %                     continue;
    %                 end
    %                 for state_idx = 1:n_states
    %                     for align_idx = 1:1
    %                         data_mat = J_data{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx};
    %                         error_mat = J_data_err{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx};
    %                         nan_filter = nan_filter & ~isnan(data_mat) & ~isnan(error_mat);
    %                         significant_filter = significant_filter | abs(data_mat) > error_mat; % criterion
    %                     end
    %                 end
    %             end
    %             filter{session_idx} = nan_filter & significant_filter;
    %         end

    %         % calculate mean and error
    %         for prepost_idx = 1:2
    %             if prepost_idx == 2 && (i == 2 || j == 2) % skip post sessions for Thalamus
    %                 continue;
    %             end
    %             for state_idx = 1:n_states
    %                 for align_idx = 1:1
    %                     data = [];
    %                     error = [];
    %                     for session_idx = 1:n_session
    %                         filter_session = filter{session_idx};
    %                         data_session = J_data{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx}(filter_session);
    %                         error_session = J_data_err{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx}(filter_session);
    %                         data = [data; data_session(:)];
    %                         error = [error; error_session(:)];
    %                     end

    %                     % nan_filter = ~isnan(data) & ~isnan(error);
    %                     % significant_filter = abs(data) > 2*error;
    %                     % % filter = nan_filter & significant_filter;
    %                     % filter = nan_filter;
    %                     data = abs(data);

    %                     % Calculate the arithmetic mean
    %                     J_mean = mean(data);

    %                     % Variance of the data points (sample standard deviation)
    %                     data_variance = var(data, 1);
    %                     sem_data = sqrt(data_variance / length(data)); % Standard error of the mean

    %                     % Propagated error from individual measurement errors
    %                     propagated_error = sqrt(sum(error.^2) / length(data)^2);

    %                     % Combine the two errors
    %                     combined_error = sqrt(sem_data^2 + propagated_error^2);

    %                     means(prepost_idx, state_idx, align_idx) = J_mean;
    %                     errors(prepost_idx, state_idx, align_idx) = combined_error;
    %                 end
    %             end
    %         end

    %         bar_width = 0.24;

    %         colors = [0,0.4470,0.7410; 0.8500,0.3250,0.0980; 0.9290,0.6940,0.1250];
    %         % histogram
    %         for align_idx = 1:1
    %             for state_idx = 1:n_states
    %                 % shift = (state_idx - (n_states+1)/2) * bar_width;
    %                 % shift = -0.5 + 0.09 + 0.29*(state_idx-1);
    %                 shift = -0.5 + 0.09 + 0.29*(state_idx-1) + 0.08*(align_idx-1);
    %                 hold on;
    %                 bar((1:2)+shift, means(:, state_idx, align_idx), bar_width, '', "FaceColor", colors(state_idx,:));
    %                 hold off;
    %             end
    %         end
    %         % errorbar
    %         for state_idx = 1:n_states
    %             for align_idx = 1:1
    %                 % shift = (state_idx - (n_states+1)/2) * bar_width;
    %                 shift = -0.5 + 0.09 + 0.29*(state_idx-1) + 0.08*(align_idx-1);
    %                 hold on;
    %                 errorbar((1:2)+shift, means(:, state_idx, align_idx), errors(:, state_idx), 'LineStyle', 'none', "Color", [0 0 0],"CapSize",1);
    %                 hold off;
    %             end
    %         end
    %         % % t-test
    %         % if i ~= 2 || j ~= 2
    %         %     for state_idx = 1:n_states
    %         %         for align_idx = 1:1
    %         %             data_pre = J_data{i, j, kernel_idx, state_idx, align_idx, :, 1};
    %         %             data_post = J_data{i, j, kernel_idx, state_idx, align_idx, :, 2};

    %         %             pos = 0.5 + 0.29*(state_idx-2);
    %         %             hold on;
    %         %             if p < 0.001
    %         %                 text(state_idx, 0.72+state_idx*0.05, '***', 'FontSize', 20, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    %         %                 plot([pos-0.5, pos+0.5], [0.7+state_idx*0.05, 0.7+state_idx*0.05], 'k', 'LineWidth', 1);
    %         %             elseif p < 0.01
    %         %                 text(state_idx, 0.72+state_idx*0.05, '**', 'FontSize', 20, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    %         %                 plot([pos-0.5, pos+0.5], [0.7+state_idx*0.05, 0.7+state_idx*0.05], 'k', 'LineWidth', 1);
    %         %             elseif p < 0.05
    %         %                 text(state_idx, 0.72+state_idx*0.05, '*', 'FontSize', 20, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    %         %                 plot([pos-0.5, pos+0.5], [0.7+state_idx*0.05, 0.7+state_idx*0.05], 'k', 'LineWidth', 1);
    %         %             end
    %         %             hold off;
    %         %         end
    %         %     end
    %         % end

    %         title([area_names{j}, ' to ', area_names{i}]);
    %         xticks(1:2);
    %         xticklabels({'Pre', 'Post'});
    %         ylabel('Average J');
    %         if i==2 && j==2
    %             % state_legends = {'Decision', '', 'Info', '',...
    %             %     'InfoAnti', '', 'InfoResp', '',...
    %             %     'RestClose', '', 'RestOpen', ''};
    %             % state_legends = {'before choice', '', 'choice to reward', '',...
    %             %     'choice to infocue', '', 'infocue to reward', '',...
    %             %     'eye close', '', 'eye open', ''};
    %             % state_legends = {'Offer1', 'Offer2', 'Decision', 'InfoAnti', 'InfoResp', 'Reward', 'RandomA', 'RandomB'};
    %             % state_legends = {'Offer1', '', 'Offer2', '', 'before choice', '', 'choice to info cue', '',...
    %             %  'after info cue', '', 'after reward', '', 'RandomA', '', 'RandomB', ''};
    %             % state_legends = {'Offer1', '', 'Offer2', '', 'before choice', '', 'Random short', '',...
    %             %  'Random Long', '', 'RandomA', '', 'RandomB', '', 'eyes open', '', 'eyes closed', ''};
    %             % state_legends = {'Random 1 second from each trial', '', 'Random 2 seconds from each trial', '', 'Random 5 seconds from each trial', ''};
    %             state_legends = {'Task', 'Eyes open', 'Eyes closed'};
    %             legend(state_legends);
    %         end
    %         % ylim(ylim_all_kernel{kernel_idx});
    %         ylim(ylim_all_kernel{1});
    %     end
    % end
    % title(tiles, [session_type,', Mean abs J, Kernel ', num2str(kernel_idx)]);

    % fig_folder = [root_path, 'figures/GLM/', session_type];
    % check_path(fig_folder);
    % saveas(f, [fig_folder, '/J_histograms_kernel_err_aligned_', num2str(kernel_idx), '.png']);

    %% bar plot of significant connection counts

    count_grid = zeros(3, 3, n_states, 3, n_session); % (pre:neg/none/pos, post: neg/none/pos, state, align)
    for state_idx = 1:n_states
        for session_idx = 1:n_session
            f = figure("Visible", "off","Position",[0, 0, 900, 900]);
            tiles = tiledlayout(2, 2);

            for i = 1:n_area
                for j = 1:n_area
                    if i==2||j==2
                        continue
                    end
                    fprintf('Kernel%d, i=%d, j=%d\n', kernel_idx, i, j);
                    nexttile;


                    % calculate filter: not nan in all states and significant in at least one state

                    nan_filter = true(size(J_data{i, j, kernel_idx, 1, 1, session_idx, 1}));
                    for prepost_idx = 1:2
                        if prepost_idx == 2 && (i == 2 || j == 2) % skip post sessions for Thalamus
                            continue;
                        end
                        % for align_idx = 1:1
                        %     data_mat = J_data{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx};
                        %     error_mat = J_data_err{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx};
                        %     nan_filter = nan_filter & ~isnan(data_mat) & ~isnan(error_mat);
                        % end
                        for align_idx = 1:1
                            for state_idx_nan = 1:n_states
                                data_mat = J_data{i, j, kernel_idx, state_idx_nan, align_idx, session_idx, prepost_idx};
                                error_mat = J_data_err{i, j, kernel_idx, state_idx_nan, align_idx, session_idx, prepost_idx};
                                nan_filter = nan_filter & ~isnan(data_mat) & ~isnan(error_mat);
                            end
                        end
                    end
                    filter = nan_filter;

                    % calculate grid
                    for align_idx = 1:1
                        pre_data = [];
                        pre_error = [];
                        post_data = [];
                        post_error = [];

                        filter_session = filter;
                        data_session = J_data{i, j, kernel_idx, state_idx, align_idx, session_idx, 1}(filter_session);
                        error_session = J_data_err{i, j, kernel_idx, state_idx, align_idx, session_idx, 1}(filter_session);
                        pre_data = [pre_data; data_session(:)];
                        pre_error = [pre_error; error_session(:)];
                        data_session = J_data{i, j, kernel_idx, state_idx, align_idx, session_idx, 2}(filter_session);
                        error_session = J_data_err{i, j, kernel_idx, state_idx, align_idx, session_idx, 2}(filter_session);
                        post_data = [post_data; data_session(:)];
                        post_error = [post_error; error_session(:)];

                        % grid count
                        for grid_i = 1:3
                            for grid_j = 1:3
                                if grid_i == 1
                                    filter_grid_i = pre_data < -filter_threshold*pre_error;
                                elseif grid_i == 2
                                    filter_grid_i = abs(pre_data) <= filter_threshold*pre_error;
                                else
                                    filter_grid_i = pre_data > filter_threshold*pre_error;
                                end
                                if grid_j == 1
                                    filter_grid_j = post_data < -filter_threshold*post_error;
                                elseif grid_j == 2
                                    filter_grid_j = abs(post_data) <= filter_threshold*post_error;
                                else
                                    filter_grid_j = post_data > filter_threshold*post_error;
                                end
                                count_grid(grid_i, grid_j, state_idx, align_idx, session_idx) = sum(filter_grid_i & filter_grid_j);
                            end
                        end
                    end

                    % plot the grid
                    data_grid = squeeze(sum(count_grid(:, :, state_idx, 1, :), 5));
                    imagesc(data_grid.');
                    % colormap('jet');
                    % colormap('gray');
                    % colormap('cool');
                    colorbar;
                    xticks(1:3);
                    xticklabels({'Negative', 'N.S.', 'Positive'});
                    yticks(1:3);
                    yticklabels({'Negative', 'N.S.', 'Positive'});
                    xlabel('Pre');
                    ylabel('Post');
                    set(gca, 'YDir', 'normal');
                    axis square;
                    xlim([0.5,3.5]);
                    ylim([0.5,3.5]);

                    for grid_i = 1:3
                        for grid_j = 1:3
                            if grid_i==2&&grid_j==2
                                text(grid_i, grid_j, num2str(data_grid(grid_i, grid_j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','Color','black');
                            else
                                text(grid_i, grid_j, num2str(data_grid(grid_i, grid_j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','Color','white');
                            end
                        end
                    end

                    title([area_names{j}, ' to ', area_names{i}]);
                    % xticks(1:2);
                    % ylim(ylim_all_kernel{kernel_idx});
                    % ylim([0 3500]);
                end
            end
                
            title(tiles, [session_type,', significant count grid, Kernel ', num2str(kernel_idx), ', ', states{state_idx}, ' Session', num2str(session_idx)]);

            fig_folder = [root_path, 'figures/GLM/', session_type];
            check_path(fig_folder);
            saveas(f, [fig_folder, '/J_histograms_kernel_count_grid_', num2str(kernel_idx), '_', states{state_idx}, '_', num2str(session_idx), '.png']);
        end
    end
end
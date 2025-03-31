% Use last aligned time bins
% calculate relative change of J count

% session_type = 'Muscimol';
session_type = 'Saline';
session_types = {'Muscimol', 'Saline'};
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

all_histogram_data = zeros(3, 3, 3, 2, 2, n_states, 3); % (i, j, kernel, pre_post, pos/neg, state, align)

n_area = 3;
ylim_all_kernel = {[-0.05, 0.9], [-0.15, 1], [-0.01, 0.12]};
for kernel_idx = 1:n_conn_kernel
    %% bar plot of significant connection counts
    f = figure("Visible", "off","Position",[0, 0, 900, 1600]);
    tiles = tiledlayout(2, 2);
    for j = 1:n_area
        for i = 1:n_area
            if (i == 2 || j == 2) % skip Thalamus
                continue;
            end
            fprintf('Kernel%d, i=%d, j=%d\n', kernel_idx, i, j);

            nexttile;

            % barplot
            count_pos = zeros(2, n_states, 3);
            count_neg = zeros(2, n_states, 3);

            % calculate filter: not nan in all states and significant in at least one state
            filter = cell(n_session, 1);
            area_size = 0; % total number of connections

            for session_idx = 1:n_session
                nan_filter = true(size(J_data{i, j, kernel_idx, 1, 1, session_idx, 1}));
                for prepost_idx = 1:2
                    if prepost_idx == 2 && (i == 2 || j == 2) % skip post sessions for Thalamus
                        continue;
                    end
                    for state_idx = 1:n_states
                        for align_idx = 1:1
                            data_mat = J_data{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx};
                            error_mat = J_data_err{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx};
                            nan_filter = nan_filter & ~isnan(data_mat) & ~isnan(error_mat);
                        end
                    end
                end
                filter{session_idx} = nan_filter;
                area_size = area_size + sum(nan_filter(:));
            end

            % calculate mean and error
            for prepost_idx = 1:2
                for state_idx = 1:n_states
                    for align_idx = 1:1
                        data = [];
                        error = [];
                        for session_idx = 1:n_session
                            filter_session = filter{session_idx};
                            data_session = J_data{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx}(filter_session);
                            error_session = J_data_err{i, j, kernel_idx, state_idx, align_idx, session_idx, prepost_idx}(filter_session);
                            data = [data; data_session(:)];
                            error = [error; error_session(:)];
                        end
                        
                        count_pos(prepost_idx, state_idx, align_idx) = sum(data > filter_threshold*error);
                        count_neg(prepost_idx, state_idx, align_idx) = sum(data < -filter_threshold*error);
                    end
                end
            end
            relative_diff_pos = (count_pos(2, :, :) - count_pos(1, :, :)) ./ count_pos(1, :, :);
            relative_diff_neg = (count_neg(2, :, :) - count_neg(1, :, :)) ./ count_neg(1, :, :);

            % all_histogram_data(i, j, kernel_idx, :, :, :) = [relative_diff_pos; relative_diff_neg];
            all_histogram_data(i, j, kernel_idx, 1, :, :, :) = count_pos/area_size;
            all_histogram_data(i, j, kernel_idx, 2, :, :, :) = count_neg/area_size;

            all_diff = [];
            all_x = [];
            bar_colors = [];

            bar_width = 0.8;

            colors = [0,0.4470,0.7410; 0.8500,0.3250,0.0980; 0.9290,0.6940,0.1250];
            lighter_colors = 1-(1-colors)*0.25;

            for state_idx = 1:n_states
                for align_idx = 1:1
                    all_diff = [all_diff, relative_diff_pos(1, state_idx, align_idx), relative_diff_neg(1, state_idx, align_idx)];
                    all_x = [all_x, state_idx*1-0.2, state_idx*1+0.2];
                    bar_colors = [bar_colors; colors(state_idx, :); lighter_colors(state_idx, :)];
                end
            end

            b = bar(all_x, all_diff, bar_width, '', "FaceColor", 'flat');
            % set barcolors
            for k = 1:2*n_states
                b.CData(k, :) = bar_colors(k, :);
            end

            title([area_names{j}, ' to ', area_names{i}]);
            xticks(1:3);
            xticklabels({'Task', 'Eye open', 'Eye close'});
            ylabel('relative change');
            ylim([-1, 15]);
        end
    end
    title(tiles, [session_type,', significant count, Kernel ', num2str(kernel_idx)]);

    fig_folder = [root_path, 'figures/GLM/', session_type];
    check_path(fig_folder);
    saveas(f, [fig_folder, '/J_histograms_kernel_rel_aligned_', num2str(kernel_idx), '.png']);
end

data_folder = [root_path, 'GLM_data/histograms/'];
check_path(data_folder);
save([data_folder, 'J_histograms_', session_type, '.mat'], 'all_histogram_data');
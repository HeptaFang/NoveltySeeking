% Require: J_cell and J_cell_err loaded
% J_data: (area i, area j, kernel, state, align, session, prepost), each cell is a n_area_i x n_area_j matrix
% J_data_err: (area i, area j, kernel, state, align, session, prepost), each cell is a n_area_i x n_area_j matrix
% J_cell: (area i, in/out/mean, kernel, state, align, session, prepost), each cell is a n_area_i x N matrix
% J_cell_err: (area i, in/out/mean, kernel, state, align, session, prepost), each cell is a n_area_i x N matrix
% firing_rate: (area, state, align, session, prepost), each cell is a n_area_i x 1 vector

% significance_all: (area, criterian), each cell is a n_neurons x 1 logical vector
% 1: pre, 2: post, 3: diff, 4: new positive, 5: new negative

% session_type = 'Muscimol';
% session_type = 'Saline';

root_path = '../';
kernel = 'DeltaPure';
reg = 'L2=2';
epoch = '2500';
area_names = {'ACC', 'Thalamus', 'VLPFC'};
states = {'Task', 'RestOpen', 'RestClose'};
aligns = {'AlignLast'};
prepost_all = {'Pre', 'Post', 'Pre_full'}; % pre cortex, post cortex, pre full
n_states = length(states);
n_aligns = length(aligns);
n_prepost = length(prepost_all);
if strcmp(session_type, 'Muscimol')
    n_session = 10;
else
    n_session = 5;
end
n_conn_kernel = 3;
n_areas = 3;

colors = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250];
lighter_colors = colors + 0.5 * (1 - colors);

property_names = {'Average J', 'Firing Rate'};
n_property = length(property_names);

% temp storage
fr_means = zeros(3, 3, 3, 2, 2); % area, state, kernel, pre/post, low/high
fr_err = zeros(3, 3, 3, 2, 2);

for align_idx = 1:n_aligns
    align = aligns{align_idx};
    for kernel_idx = 1:n_conn_kernel
        for state_idx = 1:n_states
            state = states{state_idx};
            for property_idx = 1:1 % variables to plot: 1=average J, 2=firing rate
                property_name = property_names{property_idx};
                property_data = cell(n_areas, 2, 2); % (area, low/high, prepost)
                X_all = cell(n_areas, 2); % (area, prepost) 
                Y_all = cell(n_areas, 2); % (area, prepost)
                significance_all = cell(n_areas, 5); % (area, pre/post/diff)
                fprintf('Plotting Kernel %d, %s, %s, %s\n', kernel_idx, state, align, property_name);

                for session_idx = 1:(n_session+1)
                    f = figure("Position", [100, 100, 1200, 1200], "Visible", "off");
                    t = tiledlayout(2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
                    plot_data = zeros(2, 2, 2);
                    if session_idx < n_session + 1 % Plot a single session
                        for area_i = 1:3
                            if area_i == 2 
                                continue;
                            end

                            J_pre = J_cell{area_i, 4, kernel_idx, state_idx, align_idx, session_idx, 1};
                            J_post = J_cell{area_i, 4, kernel_idx, state_idx, align_idx, session_idx, 2};
                            J_pre_err = J_cell_err{area_i, 4, kernel_idx, state_idx, align_idx, session_idx, 1};
                            J_post_err = J_cell_err{area_i, 4, kernel_idx, state_idx, align_idx, session_idx, 2};
                            firing_rate_pre = firing_rate{area_i, state_idx, align_idx, session_idx, 1};
                            firing_rate_post = firing_rate{area_i, state_idx, align_idx, session_idx, 2};
                            ave_J_pre = mean(J_pre, 2, "omitmissing");
                            ave_J_pre_err = sqrt(sum(J_pre_err.^2, 2)) / size(J_pre, 2);
                            ave_J_post = mean(J_post, 2, "omitmissing");
                            ave_J_post_err = sqrt(sum(J_post_err.^2, 2)) / size(J_post, 2);
                            ave_delta_J = ave_J_post - ave_J_pre;
                            ave_delta_J_err = sqrt(ave_J_pre_err.^2 + ave_J_post_err.^2);

                            significance_filter_pre = any(J_pre < -J_pre_err, 2);
                            significance_filter_post = any(J_post < -J_post_err, 2);
                            significance_filter_new_pos = any((J_post > J_post_err) & (J_pre <= J_pre_err), 2);
                            significance_filter_new_neg = any((J_post < -J_post_err) & (J_pre >= -J_pre_err), 2);
                            % significance_filter_post = abs(ave_J_post) > ave_J_post_err;
                            significance_filter_diff = abs(ave_delta_J) > ave_delta_J_err;
                            significance_all{area_i, 1} = [significance_all{area_i, 1}; significance_filter_pre];
                            significance_all{area_i, 2} = [significance_all{area_i, 2}; significance_filter_post];
                            significance_all{area_i, 3} = [significance_all{area_i, 3}; significance_filter_diff];
                            significance_all{area_i, 4} = [significance_all{area_i, 4}; significance_filter_new_pos];
                            significance_all{area_i, 5} = [significance_all{area_i, 5}; significance_filter_new_neg];

                            % Pre vs Delta J
                            nexttile;
                            X = firing_rate_pre;
                            % X = firing_rate_post;
                            Y = ave_delta_J;
                            X_all{area_i, 1} = [X_all{area_i, 1}; X];
                            Y_all{area_i, 1} = [Y_all{area_i, 1}; Y];
                            
                            hold on
                            scatter(X(significance_filter_diff), Y(significance_filter_diff), 30, colors(area_i, :), 'filled');
                            scatter(X(~significance_filter_diff), Y(~significance_filter_diff), 30, lighter_colors(area_i, :), 'filled');
                            if ~isempty(X)
                                plot([min(X), max(X)], [0, 0], 'k--');
                            end
                            hold off

                            xlabel('Pre Firing Rate (Hz)');
                            ylabel('Average \Delta J');
                            title([area_names{area_i}, ' ', 'Pre firing rate vs \Delta J']);

                            % Post vs Delta J
                            nexttile;
                            X = firing_rate_post;
                            % X = firing_rate_pre;
                            Y = ave_delta_J;
                            
                            X_all{area_i, 2} = [X_all{area_i, 2}; X];
                            Y_all{area_i, 2} = [Y_all{area_i, 2}; Y];

                            hold on
                            scatter(X(significance_filter_diff), Y(significance_filter_diff), 30, colors(area_i, :), 'filled');
                            scatter(X(~significance_filter_diff), Y(~significance_filter_diff), 30, lighter_colors(area_i, :), 'filled');
                            if ~isempty(X)
                                plot([min(X), max(X)], [0, 0], 'k--');
                            end
                            hold off

                            xlabel('Post Firing Rate (Hz)');
                            ylabel('Average \Delta J');
                            title([area_names{area_i}, ' ', 'Post firing rate vs \Delta J']);


                        end
                        % global title
                        sgtitle([session_type, ' Kernel ', num2str(kernel_idx), ' ', state, ' ', align, ' ', num2str(session_idx)]);
                        % saveas(f, [root_path, 'figures/deltaJ/sessions/', session_type, '_', num2str(kernel_idx), '_', state, '_', align, '_', num2str(session_idx), '.png']);
                    
                    else % Plot all sessions
                        for area_i = 1:3
                            if area_i == 2 
                                continue;
                            end
                            % Pre vs Delta J
                            nexttile;
                            X = X_all{area_i, 1};
                            Y = Y_all{area_i, 1};
                            significance_filter_pre = logical(significance_all{area_i, 1});
                            significance_filter_post = logical(significance_all{area_i, 2});
                            significance_filter_diff = logical(significance_all{area_i, 3});
                            significance_filter_new_pos = logical(significance_all{area_i, 4});
                            significance_filter_new_neg = logical(significance_all{area_i, 5});

                            hold on
                            scatter(X, Y, 30, colors(area_i, :), 'filled');
                            % scatter(X(significance_filter_diff), Y(significance_filter_diff), 30, colors(area_i, :), 'o');
                            % scatter(X(~significance_filter_diff), Y(~significance_filter_diff), 30, lighter_colors(area_i, :), 'o');
                            % % scatter(X(significance_filter_pre), Y(significance_filter_pre), 30, colors(area_i, :), 'x');
                            % % scatter(X(significance_filter_post), Y(significance_filter_post), 30, colors(area_i, :), '+');
                            % scatter(X(significance_filter_new_pos), Y(significance_filter_new_pos), 30, colors(area_i, :), 'x');
                            % scatter(X(significance_filter_new_neg), Y(significance_filter_new_neg), 30, colors(area_i, :), '+');

                            % statistics

                            ave_fr_neg = mean(X(significance_filter_new_neg));
                            err_fr_neg = std(X(significance_filter_new_neg))/sqrt(sum(significance_filter_new_neg));
                            ave_fr_pos = mean(X(significance_filter_new_pos));
                            err_fr_pos = std(X(significance_filter_new_pos))/sqrt(sum(significance_filter_new_pos));
                            ave_fr_ns = mean(X(~(significance_filter_new_neg | significance_filter_new_pos)));
                            err_fr_ns = std(X(~(significance_filter_new_neg | significance_filter_new_pos))) / sqrt(sum(~(significance_filter_new_neg | significance_filter_new_pos)));
                            
                            fprintf('Area %s, New neg: %.2f +- %.2f, New pos: %.2f +- %.2f, Non-significant: %.2f +- %.2f\n', area_names{area_i}, ave_fr_neg, err_fr_neg, ave_fr_pos, err_fr_pos, ave_fr_ns, err_fr_ns);

                            % higher 50% of delta J
                            filter_high = Y > median(Y(significance_filter_diff));

                            ave_fr_high = mean(X(filter_high));
                            err_fr_high = std(X(filter_high))/sqrt(sum(filter_high));
                            ave_fr_low = mean(X(~filter_high));
                            err_fr_low = std(X(~filter_high))/sqrt(sum(~filter_high));

                            fprintf('Area %s, Pre, High delta J: %.2f +- %.2f, Low delta J: %.2f +- %.2f\n', area_names{area_i}, ave_fr_high, err_fr_high, ave_fr_low, err_fr_low);

                            fr_means(area_i, state_idx, kernel_idx, 1, 1) = ave_fr_low;
                            fr_means(area_i, state_idx, kernel_idx, 1, 2) = ave_fr_high;
                            fr_err(area_i, state_idx, kernel_idx, 1, 1) = err_fr_low;
                            fr_err(area_i, state_idx, kernel_idx, 1, 2) = err_fr_high;
                            
                            % ave_fr_diff = mean(X(significance_filter_diff));
                            % err_fr_diff = std(X(significance_filter_diff))/sqrt(sum(significance_filter_diff));
                            % ave_fr_ns = mean(X(~significance_filter_diff));
                            % err_fr_ns = std(X(~significance_filter_diff))/sqrt(sum(~significance_filter_diff));

                            % fprintf('Area %s, Diff Significant: %.2f +- %.2f, Non-significant: %.2f +- %.2f\n', area_names{area_i}, ave_fr_diff, err_fr_diff, ave_fr_ns, err_fr_ns);

                            if ~isempty(X)
                                plot([min(X), max(X)], [0, 0], 'k--');
                            end
                            hold off

                            xlabel('Pre Firing Rate (Hz)');
                            ylabel('Average \Delta J');
                            % legend('Significant \Delta J', 'Non-significant \Delta J', 'Significant Pre', 'Significant Post');
                            % legend('Significant \Delta J', 'Non-significant \Delta J', 'New pos J', 'New neg J')
                            title([area_names{area_i}, ' ', 'Pre firing rate vs \Delta J']);

                            %% Post vs Delta J
                            nexttile;
                            X = X_all{area_i, 2};
                            Y = Y_all{area_i, 2};

                            hold on
                            scatter(X, Y, 30, colors(area_i, :), 'filled');

                            % scatter(X(significance_filter_diff), Y(significance_filter_diff), 30, colors(area_i, :), 'o');
                            % scatter(X(~significance_filter_diff), Y(~significance_filter_diff), 30, lighter_colors(area_i, :), 'o');
                            % % scatter(X(significance_filter_pre), Y(significance_filter_pre), 30, colors(area_i, :), 'x');
                            % % scatter(X(significance_filter_post), Y(significance_filter_post), 30, colors(area_i, :), '+');
                            % scatter(X(significance_filter_new_pos), Y(significance_filter_new_pos), 30, colors(area_i, :), 'x');
                            % scatter(X(significance_filter_new_neg), Y(significance_filter_new_neg), 30, colors(area_i, :), '+');
                            if ~isempty(X)
                                plot([min(X), max(X)], [0, 0], 'k--');
                            end
                            hold off

                            xlabel('Post Firing Rate (Hz)');
                            ylabel('Average \Delta J');
                            % legend('Significant \Delta J', 'Non-significant \Delta J', 'Significant Pre', 'Significant Post');
                            % legend('Significant \Delta J', 'Non-significant \Delta J', 'New pos J', 'New neg J')
                            title([area_names{area_i}, ' ', 'Post firing rate vs \Delta J']);
                            
                            % statistics
                            ave_fr_high = mean(X(filter_high));
                            err_fr_high = std(X(filter_high))/sqrt(sum(filter_high));
                            ave_fr_low = mean(X(~filter_high));
                            err_fr_low = std(X(~filter_high))/sqrt(sum(~filter_high));

                            fprintf('Area %s, Post, High delta J: %.2f +- %.2f, Low delta J: %.2f +- %.2f\n', area_names{area_i}, ave_fr_high, err_fr_high, ave_fr_low, err_fr_low);
                            fr_means(area_i, state_idx, kernel_idx, 2, 1) = ave_fr_low;
                            fr_means(area_i, state_idx, kernel_idx, 2, 2) = ave_fr_high;
                            fr_err(area_i, state_idx, kernel_idx, 2, 1) = err_fr_low;
                            fr_err(area_i, state_idx, kernel_idx, 2, 2) = err_fr_high;

                        end
                        % global title
                        sgtitle([session_type, ' Kernel ', num2str(kernel_idx), ' ', state, ' ', align, ' All']);
                        saveas(f, [root_path, 'figures/deltaJ/', session_type, '_', num2str(kernel_idx), '_', state, '_', align, '_all.png']);
                    end
                end
            end
        end
    end
end
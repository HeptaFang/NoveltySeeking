% Use last aligned time bins
% calculate relative change of J count
% OMG I must do refactoring

session_types = {'Saline', 'Muscimol'};
root_path = '../';
kernel = 'DeltaPure';
reg = 'L2=2';
epoch = '2500';
area_names = {'ACC', 'Thalamus', 'VLPFC'};
filter_threshold = 1;
n_conn_kernel = 3;

plot_mode = 'relative'; % 'prepost' or 'relative'

% load data
% states = {'Offer1', 'Offer2', 'Decision', 'InfoAnti', 'InfoResp', 'Reward', 'RandomA', 'RandomB'};
% states = {'Offer1', 'Offer2', 'Decision', 'RandomShort', 'RandomLong', 'RandomA', 'RandomB', 'RestOpen', 'RestClose'};
% states = {'RandomA', 'RandomShort', 'RandomLong'};
% states = {'Simulated', 'Simulated_higher', 'RandomLong'};
states = {'Task', 'RestOpen', 'RestClose'};
aligns = {'AlignLast'};
n_states = length(states);

all_all_histogram_data = zeros(2, 3, 3, 3, 2, 2, n_states, 3); % (session_type, i, j, kernel, pre/post, pos/neg, state, align)

for session_type_idx = 1:2 
    session_type = session_types{session_type_idx};
    data_folder = [root_path, 'GLM_data/histograms/'];
    load([data_folder, 'J_histograms_', session_type, '.mat'], 'all_histogram_data');
    all_all_histogram_data(session_type_idx, :, :, :, :, :, :, :) = all_histogram_data;
end

n_area = 3;
ylim_all_kernel = {[-0.05, 0.9], [-0.15, 1], [-0.01, 0.12]};
for kernel_idx = 1:n_conn_kernel
    %% bar plot of significant connection counts
    f = figure("Visible", "off","Position",[0, 0, 900, 600]);
    tiles = tiledlayout(2, 2, "TileSpacing", "compact", "Padding", "compact");
    for j = 1:n_area
        for i = 1:n_area
            if (i == 2 || j == 2) % skip Thalamus
                continue;
            end
            fprintf('Kernel%d, i=%d, j=%d\n', kernel_idx, i, j);
            nexttile;

            % colors = [0,0.4470,0.7410; 0.8500,0.3250,0.0980; 0.9290,0.6940,0.1250];
            colors = [0.25,0.25,0.25; 0.95,0.25,0.50];
            lighter_colors = 1-(1-colors)*0.25;

            relative_diff_pos = all_all_histogram_data(session_type_idx, i, j, kernel_idx, 1, :, :);
            relative_diff_neg = all_all_histogram_data(session_type_idx, i, j, kernel_idx, 2, :, :);

            all_bars = [];
            all_x = [];
            bar_colors = [];
            bar_width = 0.8;

            for state_idx = 1:n_states
                for session_type_idx = 1:2
                    for align_idx = 1:1
                        if strcmp(plot_mode, 'prepost')
                            for prepost_idx = 1:2
                                bar_height_pos = all_all_histogram_data(session_type_idx, i, j, kernel_idx, 1, prepost_idx, state_idx, align_idx);
                                bar_height_neg = all_all_histogram_data(session_type_idx, i, j, kernel_idx, 2, prepost_idx, state_idx, align_idx);
                                bar_pos = state_idx + (session_type_idx-1.5) * 0.35 + (prepost_idx - 1.5) * 0.15;
                                all_bars = [all_bars; bar_height_pos, bar_height_neg];
                                all_x = [all_x; bar_pos];
                                % if posneg_idx == 1
                                %     bar_colors = [bar_colors; colors(session_type_idx, :)];
                                % else
                                %     bar_colors = [bar_colors; lighter_colors(session_type_idx, :)];
                                % end
                            end
                        else
                            for posneg_idx = 1:2
                                bar_height_pre = all_all_histogram_data(session_type_idx, i, j, kernel_idx, posneg_idx, 1, state_idx, align_idx);
                                bar_height_post = all_all_histogram_data(session_type_idx, i, j, kernel_idx, posneg_idx, 2, state_idx, align_idx);
                                % relative_diff = (bar_height_post - bar_height_pre) ./ bar_height_pre;
                                relative_diff = bar_height_post ./ bar_height_pre;
                                bar_pos = state_idx + (session_type_idx-1.5) * 0.35 + (posneg_idx - 1.5) * 0.15;
                                all_bars = [all_bars; relative_diff];
                                all_x = [all_x, bar_pos];
                                if posneg_idx == 1
                                    bar_colors = [bar_colors; colors(session_type_idx, :)];
                                else
                                    bar_colors = [bar_colors; lighter_colors(session_type_idx, :)];
                                end
                            end

                        end
                    end
                end
            end

            b = bar(all_x, all_bars.', bar_width, 'stacked', "FaceColor", 'flat');

            % 4 invisible auxilary bars for legend
            hold on;
            bar(1, 0, 'FaceColor', colors(1, :));
            bar(1, 0, 'FaceColor', lighter_colors(1, :));
            bar(1, 0, 'FaceColor', colors(2, :));
            bar(1, 0, 'FaceColor', lighter_colors(2, :));
            hold off;

            title([area_names{j}, ' to ', area_names{i}]);
            xticks(1:3);
            xticklabels({'Task', 'Eye open', 'Eye close'});
            if strcmp(plot_mode, 'prepost')
                ylabel('significant connection ratio');
                % set barcolors
                for state_idx = 1:n_states
                    for session_type_idx = 1:2
                        for prepost_idx = 1:2
                            k = (state_idx-1)*4 + (session_type_idx-1)*2 + prepost_idx;
                            b(1).CData(k, :) = colors(session_type_idx, :);
                            b(2).CData(k, :) = lighter_colors(session_type_idx, :);
                        end
                    end
                end
                % ylim([-1, 15]);
                xlim([0.6, 3.4]);
                ylim([0, 0.24]);
                legend({'', '', 'Saline positive', 'Saline negative', 'Muscimol positive', 'Muscimol negative'}, 'Location', 'northwest');
            else
                % ylabel('relative change');
                ylabel('post/pre');
                for state_idx = 1:n_states
                    for session_type_idx = 1:2
                        for posneg_idx = 1:2
                            k = (state_idx-1)*4 + (session_type_idx-1)*2 + posneg_idx;
                            b.CData(k, :) = bar_colors(k, :);
                        end
                    end
                end

                hold on;
                plot([0,4], [1, 1], 'k--');
                hold off;
                ylim([0, 15]);
                xlim([0.6, 3.4]);
                % ylim([0, 0.24]);
                legend({'', 'Saline positive', 'Saline negative', 'Muscimol positive', 'Muscimol negative'}, 'Location', 'north');
            end

        end
    end
    title(tiles, [session_type,', significant count, Kernel ', num2str(kernel_idx)]);

    fig_folder = [root_path, 'figures/GLM/histogram'];
    check_path(fig_folder);
    saveas(f, [fig_folder, '/J_histograms_', plot_mode, '_kernel', num2str(kernel_idx), '.png']);
end

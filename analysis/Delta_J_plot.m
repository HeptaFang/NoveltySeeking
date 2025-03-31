% Req: temp storage
% fr_means = zeros(3, 3, 3, 2, 2); % area, state, kernel, pre/post, low/high
% fr_err = zeros(3, 3, 3, 2, 2);

% session_type = 'Muscimol';
n_areas = 2;
n_states = 3;
n_conn_kernel = 3;

area_names = {'ACC', 'thalamus', 'VLPFC'};
states = {'Task', 'EyesOpen', 'EyesClose'};

for kernel_idx = 1:n_conn_kernel
    f = figure("Position", [100, 100, 1200, 800], "Visible", "off");
    t = tiledlayout(2, 3);
    for i = 1:n_areas
        for state_idx = 1:n_states
            nexttile(i*3+state_idx-3)
            if i == 2
                i_eff = 3;
            else
                i_eff = i;
            end
            hold on;
            data = fr_means(i_eff, state_idx, kernel_idx, :, :);
            % swap axes to (low/high, pre/post)
            data = permute(data, [5, 4, 1, 2, 3]);
            bar(1:2, data);

            % error bars
            err = fr_err(i_eff, state_idx, kernel_idx, :, :); 
            err = permute(err, [5, 4, 1, 2, 3]);
            for j = 1:2
                for k = 1:2
                    errorbar(j + (0.3 * (k-1.5)), data(j, k), err(j, k), 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
                end
            end

            hold off;
            % set axis limits and labels
            xlim([0.5, 2.5]);
            xticks(1:2);
            xticklabels({'Pre', 'Post'});
            max_y = max(data(:) + err(:));
            ylim([0, max_y * 1.25]);
            title([area_names{i_eff}, ' ', states{state_idx}]);
            xlabel('Session Type');
            ylabel('Firing Rate');
            legend({'Low \Delta J', 'High \Delta J'}, 'Location', 'north');
        end
    end
    % set figure title and save
    title(t, sprintf('%s Kernel %d', session_type, kernel_idx), 'FontSize', 16, 'FontWeight', 'bold');
    saveas(f, ['../figures/deltaJ/Delta_J_plot_', session_type, '_kernel_', num2str(kernel_idx), '.png']);
end
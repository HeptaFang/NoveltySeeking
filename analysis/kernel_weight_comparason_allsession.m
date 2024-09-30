function kernel_weight_comparason_allsession(all_diffs)


%%%%%% histogram of kernel weight differences

area_names = {'ACC', 'Thalamus', 'VLPFC'};
fig = figure("Visible", "off");
set(fig, 'PaperPosition', [0, 0, 8, 6]);
tiles = tiledlayout(3, 3);
for i = 1:3
    for j = 1:3
        ax = nexttile;
        
        % Plot each pair of elements as a line
        hold on;
        diff = all_diffs{i, j};
        h = histogram(diff(:), 40, "Normalization", "count", "BinLimits", [-1, 1]);
        xlim([-1, 1]);
        ylim([0, 1.3 * max(h.Values)]);
        xline(0, 'k-', 'LineWidth', 1.0);

        % Plot vertical lines for mean and median
        mean_diff = mean(diff(:));
        % median_diff = median(diff(:));
        if mean_diff > 0
            xline(mean_diff, 'r--', 'Label', 'Mean>0', 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'FontSize', 8);
        else
            xline(mean_diff, 'b--', 'Label', 'Mean<0', 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', 'FontSize', 8);
        end
        % xline(median_diff, 'b--', 'Label', 'Median', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');

        % test if the mean is significantly different from 0
        [~, p] = ttest(diff(:), 0, 'Alpha', 0.001);
        if p < 0.001
            text(0.5, 5/6, '***', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'Color', 'r');
        elseif p < 0.01
            text(0.5, 5/6, '**', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'Color', 'r');
        elseif p < 0.05
            text(0.5, 5/6, '*', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'Color', 'r');
        end

        % Set plot properties
        % title(['Area ' num2str(i) ' vs ' num2str(j)]);
        title([area_names{i} ' ' area_names{j}]);
        xlabel('Fast-Slow');
        ylabel('count');
        grid on;
    end
end

% Add a super title for the entire figure
sgtitle("All sessions");

% Save the figure
fig_path = ['../figures/GLM'];
check_path(fig_path);
fig_file = [fig_path, '/GLMkernels_hist_all_sorted.png'];
print(fig, fig_file, '-dpng', '-r200');

end

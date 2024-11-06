function plot_density_histogram()
    % Load the data
    all_diffs = load('../GLM_data/all_diffs.mat');
    all_Js = load('../GLM_data/all_Js.mat');
    
    % Define figure and subplots (3x3 grid)
    figure;
    area_names = {'ACC', 'Thalamus', 'VLPFC'};
    
    max_log_density = -inf;  % To keep track of max log density
    
    for i = 1:3
        for j = 1:3
            fast_Js = all_Js.all_Js{i, j, 1}(:);  % Flatten
            slow_Js = all_Js.all_Js{i, j, 2}(:);  % Flatten
            
            % Create a 2D histogram
            [hist_counts, xedges, yedges] = histcounts2(fast_Js, slow_Js, 100, 'XBinLimits', [-2, 2]);
            hist_log = log1p(hist_counts);  % log1p to avoid log(0)
            max_log_density = max(max_log_density, max(hist_log(:)));  % Track max log density
            
            % Plot the 2D histogram
            subplot(3, 3, (i-1)*3 + j);
            imagesc(xedges, yedges, hist_log');
            hold on;
            plot([-2, 2], [-2, 2], 'k--', 'LineWidth', 1);  % Diagonal line
            
            axis equal;
            xlim([-2 2]);
            ylim([-2 2]);
            set(gca, 'XTick', [-2 0 2], 'YTick', [-2 0 2]);
            
            % Add gray lines at zero
            line([0, 0], [-2, 2], 'Color', 'black', 'LineWidth', 1);
            line([-2, 2], [0, 0], 'Color', 'black', 'LineWidth', 1);
            
            % Overlay histograms along x and y axes
            ax_hist_x = axes('Position', get(gca, 'Position') + [0, -0.1, 0, 0]); 
            histogram(fast_Js, 50, 'BinLimits', [-2, 2], 'FaceColor', 'k');
            set(ax_hist_x, 'XTick', [], 'YTick', []);
            ax_hist_y = axes('Position', get(gca, 'Position') + [-0.1, 0, 0, 0]);
            histogram(slow_Js, 50, 'BinLimits', [-2, 2], 'Orientation', 'horizontal', 'FaceColor', 'k');
            set(ax_hist_y, 'XTick', [], 'YTick', []);
            
            % Labeling
            xlabel('J\_fast');
            ylabel('J\_slow');
            title([area_names{j} ' to ' area_names{i}]);
        end
    end
    
    % Colorbar setup
    c = colorbar('Location', 'eastoutside');
    clim([0, max_log_density]);
    c.Label.String = 'Density';
    
    % Add global title
    sgtitle('Density plots of J\_fast vs J\_slow (abs J)');
    
    % Save figure
    saveas(gcf, 'figures/kernel_comp_density_abs.png');
end

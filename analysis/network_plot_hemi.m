function network_plot_hemi(ax, mat, err, borders, node_colors, omitThr, area_mode, dist_mode)
    % Plot three areas in a circle
    % ACC: 1-borders(1), Thalamus: borders(1)-borders(2), VLPFC: borders(2)-N

    % default values
    if nargin < 7
        area_mode = 'full'; % Included areas: full or cortex
    end
    if nargin < 8
        dist_mode = 'full'; % Distribution: area or cell
    end

    borders = borders(2:end);

    % N = size(mat, 1);
    N_ACC = borders(1) - 1;
    N_thalamus = borders(2) - borders(1);
    N_VLPFC = borders(3) - borders(2);

    if strcmp(area_mode, 'full')
        N = N_ACC + N_thalamus + N_VLPFC;

        theta_ACC = linspace((1/2)*pi, (7/6)*pi, N_ACC);
        theta_thalamus = linspace((7/6)*pi, (11/6)*pi, N_thalamus);
        theta_VLPFC = linspace((11/6)*pi, (15/6)*pi, N_VLPFC);

        % ACC: top left, VLPFC: top right, Thalamus: lower hemisphere
        x_ACC = cos(theta_ACC) + 0.5*cos((5/6)*pi);
        y_ACC = sin(theta_ACC) + 0.5*sin((5/6)*pi);
        x_thalamus = cos(theta_thalamus) + 0.5*cos((9/6)*pi);
        y_thalamus = sin(theta_thalamus) + 0.5*sin((9/6)*pi);
        x_VLPFC = cos(theta_VLPFC) + 0.5*cos((13/6)*pi);
        y_VLPFC = sin(theta_VLPFC) + 0.5*sin((13/6)*pi);
        x = [x_ACC, x_thalamus, x_VLPFC];
        y = [y_ACC, y_thalamus, y_VLPFC];
        hold(ax, "on");
        
        for i = 1:N
            for j = 1:N
                if ((i >= borders(1) && i < borders(2)) || (j >= borders(1) && j < borders(2))) && omitThr
                    continue;
                end
                if mat(i, j) > err(i, j)
                    % red solid line
                    quiver(ax, x(i), y(i), x(j) - x(i), y(j) - y(i), 0, 'r', 'LineStyle', '-');
                elseif mat(i, j) < -err(i, j)
                    % blue dashed line
                    quiver(ax, x(i), y(i), x(j) - x(i), y(j) - y(i), 0, 'b', 'LineStyle', '--');
                end
            end
        end
        for i = 1:N
            area_init = 0;
            if i >= borders(2)
                area_init = borders(2) - 1;
            end
            if i >= borders(3)
                area_init = borders(3) - 1;
            end
            plot(ax, x(i), y(i), 'o', 'MarkerSize', 3, 'MarkerFaceColor', node_colors(i, :), 'MarkerEdgeColor', 'k');
            text(ax, x(i)*1.1, y(i)*1.1, num2str(i-area_init), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        % for i = (border + 1):N
        %     plot(ax, x(i), y(i), 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
        %     text(ax, x(i)*1.1, y(i)*1.1, num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        % end
        hold(ax, "off");
        axis(ax, "equal");
        xlim(ax, [-1.5, 1.5]);
        ylim(ax, [-1.5, 1.5]);
        xticks(ax, []);
        yticks(ax, []);

        ax.XAxis.Visible = 'off';
        ax.YAxis.Visible = 'off';

        text(ax, -1.5, 1.5, 'ACC', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
        text(ax, 0, -1.4, 'Thalamus', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
        text(ax, 1.5, 1.5, 'VLPFC', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    else % cortex
        N = N_ACC + N_VLPFC;
        theta_ACC = linspace((4/6)*pi, (8/6)*pi, N_ACC);
        theta_VLPFC = linspace((10/6)*pi, (14/6)*pi, N_VLPFC);

        % ACC: left, VLPFC: right
        x_ACC = cos(theta_ACC);
        y_ACC = sin(theta_ACC);
        x_VLPFC = cos(theta_VLPFC);
        y_VLPFC = sin(theta_VLPFC);
        x = [x_ACC, x_VLPFC];
        y = [y_ACC, y_VLPFC];

        hold(ax, "on");
        for i = 1:N
            for j = 1:N
                if ((i >= borders(1) && i < borders(2)) || (j >= borders(1) && j < borders(2)))
                    continue;
                end
                if mat(i, j) > err(i, j)
                    % red solid line
                    quiver(ax, x(i), y(i), x(j) - x(i), y(j) - y(i), 0, 'r', 'LineStyle', '-');
                elseif mat(i, j) < -err(i, j)
                    % blue dashed line
                    quiver(ax, x(i), y(i), x(j) - x(i), y(j) - y(i), 0, 'b', 'LineStyle', '--');
                end
            end
        end

        for i = 1:N
            area_init = 0;
            if i >= borders(2)
                area_init = borders(2) - 1;
            end
            plot(ax, x(i), y(i), 'o', 'MarkerSize', 3, 'MarkerFaceColor', node_colors(i, :), 'MarkerEdgeColor', 'k');
            text(ax, x(i)*1.1, y(i)*1.1, num2str(i-area_init), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        % for i = (border + 1):N
        %     plot(ax, x(i), y(i), 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
        %     text(ax, x(i)*1.1, y(i)*1.1, num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        % end
        hold(ax, "off");
        axis(ax, "equal");
        xlim(ax, [-1, 1]);
        ylim(ax, [-1, 1]);
        xticks(ax, []);
        yticks(ax, []);
        ax.XAxis.Visible = 'off';
        ax.YAxis.Visible = 'off';
        text(ax, -1, 1, 'ACC', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
        text(ax, 1, 1, 'VLPFC', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    end
end
function [all_diffs, all_Js] = kernel_weight_comparason(dataset_name, session, kernel_name, epoch, reg, shuffle_size, sorting, idx)
% Load the original model parameters and data
model_path_ori = ['../GLM_model/', dataset_name, '/GLM_', dataset_name, '_', ...
        int2str(session), '_', kernel_name, '_0_', ...
        reg.name, '_', int2str(epoch)];
load(model_path_ori, "model_par", "PS_kernels", "conn_kernels", "n_PS_kernel", "n_conn_kernel", "kernel_len", "N");

data_path = ['../GLM_data/', dataset_name, '/GLMdata_', dataset_name, '_', ...
        int2str(session),'_', kernel_name, '_0.mat'];
load(data_path, "raster");
firing_rate = mean(raster, 2);

% Determine session border
if session<100
    session_border = session;
else
    session_border = floor(session/100);
end

% Load the borders for the specific session
load(['../GLM_data/', dataset_name,'/borders_', dataset_name, '_', ...
        int2str(session_border),'.mat'], "borders");

% For acc-thalamus-vlpfc dataset
A_T_border = borders(1);
T_P_border = borders(2);
    
par_ori = model_par;

% Load shuffled parameters for statistical testing
par_sfl = zeros([size(par_ori), shuffle_size]);
for i = 1:shuffle_size
    model_path_sfl = ['../GLM_model/', dataset_name, '/GLM_', dataset_name, '_', ...
        int2str(session), '_', kernel_name, '_', int2str(i), '_', ...
        reg.name, '_', int2str(epoch)];
    par_sfl(:, :, i) = load(model_path_sfl).model_par;
end

% Perform t-tests to find significant differences
significant = zeros(size(par_ori));
for i = 1:N
    for j = 1:(1 + n_PS_kernel + N * n_conn_kernel)
        % One-sample t-test
        h = ttest(reshape(par_sfl(i, j, :) - par_ori(i, j), [1, shuffle_size]), ...
            0, "Alpha", 0.001);
        significant(i, j) = h;
    end
end

% Calculate significant parameter differences
par_sig = par_ori;
par_sig = (par_sig - mean(par_sfl, 3)) .* significant;
save("par_sig.mat", "par_sig");

% Set limits for plotting
clim_ori = max(abs(par_ori(:, (n_PS_kernel + 2):end)), [], "all");
clim_sfl = max(abs(par_sfl(:, (n_PS_kernel + 2):end, :)), [], "all");
clim_all = max(clim_ori, clim_sfl);
clim_all = 2;

% Sorting of parameters based on defined criteria
sorting_ranges = [1, A_T_border - 0.5; A_T_border + 0.5, T_P_border - 0.5; T_P_border + 0.5, N];
criterion = 1;
sort_idx = zeros(1, N);

criterion_mat = par_sig(:, (N * (criterion - 1) + n_PS_kernel + 2):(N * criterion + n_PS_kernel + 1));
criterion_mat(isnan(criterion_mat)) = 0;

% Sort based on the given criterion
for i = 1:3
    sort_s = sorting_ranges(i, 1);
    sort_e = sorting_ranges(i, 2);

    if sorting == "count"
        s = criterion_mat ~= 0;
        inward = sum(s, 2).';
        outward = sum(s, 1);
        v = inward + outward;
        [~, sort_idx(sort_s:sort_e)] = sort(v(sort_s:sort_e), 'descend');
    elseif sorting == "abs"
        s = abs(criterion_mat);
        inward = sum(s, 2).';
        outward = sum(s, 1);
        v = inward + outward;
        [~, sort_idx(sort_s:sort_e)] = sort(v(sort_s:sort_e), 'descend');
    elseif sorting == "square"
        s = criterion_mat.^2;
        inward = sum(s, 2).';
        outward = sum(s, 1);
        v = inward + outward;
        [~, sort_idx(sort_s:sort_e)] = sort(v(sort_s:sort_e), 'descend');
    elseif sorting == "max"
        s = abs(criterion_mat);
        inward = max(s, [], 2).';
        outward = max(s, [], 1);
        v = max(inward, outward);
        [~, sort_idx(sort_s:sort_e)] = sort(v(sort_s:sort_e), 'descend');
    elseif sorting == "idx"
        [~, sort_idx(sort_s:sort_e)] = sort(idx(sort_s:sort_e), 'ascend');
    end
    sort_idx(sort_s:sort_e) = sort_idx(sort_s:sort_e) + sort_s - 1;
end

% Apply sorting to the parameters
par_ori(:, :) = par_ori(sort_idx, :);
par_sig(:, :) = par_sig(sort_idx, :);
par_sfl(:, :, :) = par_sfl(sort_idx, :, :);

for j = 1:n_conn_kernel
    shift = N * (j - 1) + n_PS_kernel + 1;
    par_ori(:, shift + 1:shift + N) = par_ori(:, shift + sort_idx);
    par_sig(:, shift + 1:shift + N) = par_sig(:, shift + sort_idx);
    par_sfl(:, shift + 1:shift + N, :) = par_sfl(:, shift + sort_idx, :);
end

% Save the sorting index
sort_path = ['../GLM_data/', dataset_name, '/sortidx_', dataset_name, '_', ...
        int2str(session), '_', kernel_name, '.mat'];
save(sort_path, "sort_idx");

area_names = {'ACC', 'Thalamus', 'VLPFC'};

kernel1 = par_ori(:, 4:3 + N);
kernel2 = par_ori(:, 4 + N:3 + 2 * N);
%% Scatter plot of kernel pairs
fig = figure("Visible", "off");
set(fig, 'PaperPosition', [0, 0, 6, 6]);
tiles = tiledlayout(3, 3);
cmap = brewermap(256,'*RdBu');

% Plot comparison between kernel1 and kernel2 across the areas
for i = 1:3
    for j = 1:3
        ax = nexttile;
        sort_s_x = sorting_ranges(i, 1);
        sort_e_x = sorting_ranges(i, 2);
        sort_s_y = sorting_ranges(j, 1);
        sort_e_y = sorting_ranges(j, 2);
    
        % Extract the relevant kernels for the area
        kernel1_area = kernel1(sort_s_x:sort_e_x, sort_s_y:sort_e_y);
        kernel2_area = kernel2(sort_s_x:sort_e_x, sort_s_y:sort_e_y);

        kernel1_area = abs(kernel1_area);
        kernel2_area = abs(kernel2_area);
        
        % Plot each pair of elements as a line
        hold on;
        scatter(kernel1_area(:), kernel2_area(:), 2, "Marker", ".", "MarkerEdgeColor", "r");
        value_lim = max([max(kernel1_area(:)), max(kernel2_area(:))]);
        plot([0, value_lim], [0, value_lim], 'k--');
        
        % Set plot properties
        % title(['Area ' num2str(i) ' vs ' num2str(j)]);
        title([area_names{i} ' ' area_names{j}]);
        if value_lim > 0
            xlim([0, value_lim]);
            ylim([0, value_lim]);
        end
        xlabel('Kernel1');
        ylabel('Kernel2');
        grid on;
    end
end

% Add a super title for the entire figure
sgtitle([dataset_name, ', session ', int2str(session), ', ', kernel_name]);

% Save the figure
fig_path = ['../figures/GLM/', dataset_name];
check_path(fig_path);
fig_file = [fig_path, '/GLMkernels_' dataset_name, '_', ...
        int2str(session), '_', kernel_name, '_', ...
        reg.name, '_', int2str(epoch), '_sorted.png'];
print(fig, fig_file, '-dpng', '-r200');

%%%%%% histogram of kernel weight differences

fig = figure("Visible", "off");
set(fig, 'PaperPosition', [0, 0, 8, 6]);
tiles = tiledlayout(3, 3);
all_diffs = cell(3, 3);
all_Js = cell(3, 3, 2);
for i = 1:3
    for j = 1:3
        ax = nexttile;
        sort_s_x = sorting_ranges(i, 1);
        sort_e_x = sorting_ranges(i, 2);
        sort_s_y = sorting_ranges(j, 1);
        sort_e_y = sorting_ranges(j, 2);
    
        % Extract the relevant kernels for the area
        kernel1_area = kernel1(sort_s_x:sort_e_x, sort_s_y:sort_e_y);
        kernel2_area = kernel2(sort_s_x:sort_e_x, sort_s_y:sort_e_y);

        kernel1_area = abs(kernel1_area);
        kernel2_area = abs(kernel2_area);
        
        % Plot each pair of elements as a line
        hold on;
        diff = kernel1_area - kernel2_area;
        all_diffs{i, j} = diff(:);
        all_Js{i, j, 1} = kernel1_area(:);
        all_Js{i, j, 2} = kernel2_area(:);
        histogram(diff(:), 40, "Normalization", "count", "BinLimits", [-1, 1]);
        xline(0, 'k-', 'LineWidth', 1.0);

        % Plot vertical lines for mean and median
        mean_diff = mean(diff(:));
        % median_diff = median(diff(:));
        if mean_diff > 0
            xline(mean_diff, 'r--', 'Label', 'Mean', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');
        else
            xline(mean_diff, 'r--', 'Label', 'Mean', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top');
        end
        % xline(median_diff, 'b--', 'Label', 'Median', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');

        % test if the mean is significantly different from 0
        [~, p] = ttest(diff(:), 0, 'Alpha', 0.001);
        if p < 0.001
            text(0.5, 0.5, '***', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'Color', 'r');
        elseif p < 0.01
            text(0.5, 0.5, '**', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'Color', 'r');
        elseif p < 0.05
            text(0.5, 0.5, '*', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'Color', 'r');
        end

        % Set plot properties
        % title(['Area ' num2str(i) ' vs ' num2str(j)]);
        title([area_names{i} ' ' area_names{j}]);
        xlabel('Fast-Slow');
        ylabel('count');
        xlim([-1, 1]);
        grid on;
    end
end

% Add a super title for the entire figure
sgtitle([dataset_name, ', session ', int2str(session), ', ', kernel_name]);

% Save the figure
fig_path = ['../figures/GLM/', dataset_name];
check_path(fig_path);
fig_file = [fig_path, '/GLMkernels_hist_' dataset_name, '_', ...
        int2str(session), '_', kernel_name, '_', ...
        reg.name, '_', int2str(epoch), '_sorted.png'];
print(fig, fig_file, '-dpng', '-r200');

end

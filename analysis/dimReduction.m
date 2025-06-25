% PCA of raster data

clear; clc; close all;
% load data
% models = {'Ctr', 'Sync', 'Mus'};
% models = {'Ctr', 'Mus'};
% session_type = 'Muscimol';
% session_type = 'Saline';
session_type = 'Simulated';
load_mode = 'GLM'; % 'raw' or 'GLM'
% load_mode = 'raw'; % 'raw' or 'GLM'

if strcmp(session_type, 'Muscimol')
    sessions = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    models = {'PreTask', 'PostTask', 'PreRestOpen', 'PostRestOpen', 'PreRestClose', 'PostRestClose', };
    linestyles = {'r-', 'r--', 'b-', 'b--', 'm-', 'm--', };

elseif strcmp(session_type, 'Saline')
    sessions = [1, 2, 3, 4, 5];
    models = {'PreTask', 'PostTask', 'PreRestOpen', 'PostRestOpen', 'PreRestClose', 'PostRestClose', };
    linestyles = {'r-', 'r--', 'b-', 'b--', 'm-', 'm--', };

elseif strcmp(session_type, 'Simulated')
    % sessions = [1];
    sessions = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    models = {'Ctr', 'Sync', 'Mus'};
    linestyles = {'r-', 'b-', 'm-',};
end

n_models = length(models);
n_sessions = length(sessions);
explain_threshold = 90;

myGauss = @(x, mu, sigma) exp(-((x-mu).^2)/(2*sigma^2));

% gaussian smoothing
% kernel_x = -500:500;
% kernel_y = myGauss(kernel_x, 0, 120);
kernel_x = -200:200;
kernel_y = myGauss(kernel_x, 0, 20);
smooth_kernel = kernel_y/sum(kernel_y); % normalize to 1
% smooth_kernel = ones(1,1000)/1000; % uniform smoothing

figure;
tiledlayout(2,5);
dim_nums = zeros(n_models, n_sessions); % model, session
dim_scores = zeros(n_models, n_sessions); % model, session
merged_dim_nums = zeros(1, n_sessions); % session
merged_dim_scores = zeros(1, n_sessions); % session

for session_idx = sessions
    if strcmp(session_type, 'Muscimol') && session_idx == 2
        continue;
    end
    % cell filter
    cell_filter = NaN;
    for model_idx = 1:n_models
        model_name = models{model_idx};
        if strcmp(session_type, 'Muscimol') || strcmp(session_type, 'Saline')
            model_name_full = [session_type, model_name, '_cortex_AlignLast'];
        elseif strcmp(session_type, 'Simulated')
            model_name_full = ['Simulated_', model_name, '_partial'];
            % model_name_full = ['Simulated_', model_name, '_partial'];
        end

        % load model
        if strcmp(load_mode, 'raw')
            file_path = ['../GLM_data/', model_name_full, '/raster_', model_name_full, '_', num2str(session_idx), '_0.mat'];
            load(file_path, 'rasters', 'n_trial', 'trial_len');
            rasters_concat = cell2mat(rasters);
        else
            file_path = ['../GLM_data/', model_name_full, '/GLMdata_', model_name_full, '_', num2str(session_idx), '_DeltaPure_0.mat'];
            load(file_path, 'raster');
            rasters_concat = raster; % load raster data
        end

        N = size(rasters_concat, 1);
        if isnan(cell_filter)
            cell_filter = logical(sum(rasters_concat, 2) > 0); % filter out silent cells
        else
            cell_filter = cell_filter & logical(sum(rasters_concat, 2) > 0); % filter out silent cells
        end
    end
    fprintf('Session %d, Number of active cells: %d/%d\n', session_idx, sum(cell_filter), N);

    explained_all = cell(n_models);
    explained_all_cumsum = cell(n_models);
    active_num = zeros(n_models, 1);
    merged_raster = cell(1, n_models);

    for model_idx = 1:n_models
        model_name = models{model_idx};
        if strcmp(session_type, 'Muscimol') || strcmp(session_type, 'Saline')
            model_name_full = [session_type, model_name, '_cortex_AlignLast'];
        elseif strcmp(session_type, 'Simulated')
            model_name_full = ['Simulated_', model_name, '_partial'];
            % model_name_full = ['Simulated_', model_name, '_partial'];
        end

        % load model
        if strcmp(load_mode, 'raw')
            file_path = ['../GLM_data/', model_name_full, '/raster_', model_name_full, '_', num2str(session_idx), '_0.mat'];
            load(file_path, 'rasters', 'n_trial', 'trial_len');
            rasters_concat = cell2mat(rasters);
        else
            file_path = ['../GLM_data/', model_name_full, '/GLMdata_', model_name_full, '_', num2str(session_idx), '_DeltaPure_0.mat'];
            load(file_path, 'raster');
            rasters_concat = raster; % load raster data
        end

        rasters_concat = rasters_concat(cell_filter, :); % filter out silent cells
        disp(size(rasters_concat)); % matrix of size (N, n_trial*trial_len)

        % smooth rasters
        rasters_concat = conv2(rasters_concat, smooth_kernel, 'valid');
        % z-score for each cell
        rasters_concat = zscore(rasters_concat, 0, 2); % z-score for each cell

        merged_raster{model_idx} = rasters_concat;

        % PCA
        [coeff,score,latent,tsquared,explained,mu] = pca(rasters_concat');
        explained_all{model_idx} = explained;
        explained_all_cumsum{model_idx} = cumsum(explained);
        active_num(model_idx) = size(rasters_concat, 1); % number of active cells
    end

    merged_raster = cell2mat(merged_raster); % merge all rasters

    % PCA
    [coeff,score,latent,tsquared,explained,mu] = pca(merged_raster');
    merged_cumsum = cumsum(explained);
    % merged_dim_num = find(merged_cumsum >= explain_threshold, 1, 'first');
    merged_dim_num = 2^calc_entropy(explained); % entropy of the number of components
    merged_dim_nums(session_idx) = merged_dim_num;
    merged_dim_scores(session_idx) = merged_dim_num/size(merged_raster, 1); % normalized by number of active cells

    nexttile;
    hold on;
    for model_idx = 1:n_models
        plot(explained_all_cumsum{model_idx}, linestyles{model_idx});
    end
    plot(merged_cumsum, 'k-');
    hold off;
    % legend(models);
    title(['Session ', num2str(session_idx)]);

    % find the number of components that explain [explain_threshold] variance
    for model_idx = 1:n_models
        explained = explained_all{model_idx};
        explained_cumsum = explained_all_cumsum{model_idx};
        % dim_num = find(explained_cumsum >= explain_threshold, 1, 'first');
        dim_num = 2^calc_entropy(explained);
        dim_nums(model_idx, session_idx) = dim_num;
        dim_scores(model_idx, session_idx) = dim_num/active_num(model_idx); % normalized by number of active cells
    end
end

if strcmp(session_type, 'Muscimol')
    % remove session 2
    dim_scores = dim_scores(:, [1, 3:end]);
    dim_nums = dim_nums(:, [1, 3:end]);
    merged_dim_scores = merged_dim_scores([1, 3:end]);
    merged_dim_nums = merged_dim_nums([1, 3:end]);
    sessions = [1, 3, 4, 5, 6, 7, 8, 9, 10];
end

% plot statistics
figure;
tiledlayout(3, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact');
nexttile;
hold on;
for model_idx = 1:n_models
    plot(dim_nums(model_idx, :), linestyles{model_idx});
end
plot(merged_dim_nums, 'k-');
hold off;
% legend(models);
legend([models, 'Merged']);
% title(['Number of components that explain ', num2str(explain_threshold), '% variance']);
title('Number of effective components');
xlabel('Session');
ylabel('Number of components');
xticklabels(sessions);
xticks(1:n_sessions);
xlim([0, 11]);


nexttile;
hold on;
for model_idx = 1:n_models
    plot(dim_scores(model_idx, :), linestyles{model_idx});
end
plot(merged_dim_scores, 'k-');
hold off;
% legend(models);
legend([models, 'Merged']);
% title(['Normalized number of components that explain ', num2str(explain_threshold), '% variance']);
title('Normalized number of effective components');
xlabel('Session');
ylabel('Normalized number of components');
xticklabels(sessions);
xticks(1:n_sessions);
xlim([0, 11]);

nexttile;
hold on;
for model_idx = 1:n_models
    % plot(dim_scores(model_idx, :)./merged_dim_scores, linestyles{model_idx});
    plot(dim_scores(model_idx, :)./dim_scores(1, :), linestyles{model_idx});
end
hold off;
% legend(models);
legend([models]);
% title(['Merge-Normalized number of components that explain ', num2str(explain_threshold), '% variance']);
title('Merge-Normalized number of effective components');
xlabel('Session');
ylabel('Merge-Normalized number of components');
xticklabels(sessions);
xticks(1:n_sessions);
xlim([0, 11]);

dim_scores = dim_scores(:, :)./dim_scores(1, :);

score_ave = mean(dim_scores, 2);
score_std = std(dim_scores, 0, 2);
score_se = score_std/sqrt(n_sessions-1); % standard error of the mean

figure;
bar(score_ave, 'FaceColor', 'flat', 'EdgeColor', 'none');
hold on;
errorbar(1:n_models, score_ave, score_se, 'k.', 'LineWidth', 2);
hold off;
% title([session_type, ', Ratio of components, ', num2str(explain_threshold), '% variance']);
title([session_type, ', Normalized effective dimension']);
xlabel('Model');
ylabel('Normalized effective dimension');
xticks(1:n_models);
xticklabels(models);
ylim([0.65, 0.85])


function e = calc_entropy(x)
    % calculate entropy of a vector
    p = x/sum(x);
    p(p == 0) = [];
    e = -sum(p.*log2(p));
end
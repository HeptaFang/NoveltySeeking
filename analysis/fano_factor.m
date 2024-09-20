% load data
controls = {'Muscimol', 'Saline'};
n_sessions = [10, 5];
% subdatasets = {'pre', 'post', 'VLPFC'};

% Create figure
figure;
hold on;

for control_idx = 1:2
    control = controls{control_idx};
    session_num = n_sessions(control_idx);

    for session_idx = 1:session_num
        %% load data
        dataset_pre = [control, 'Pre_cortex'];
        dataset_post = [control, 'Post_cortex'];
        pre_file = ['../GLM_data/', dataset_pre, '/raster_', dataset_pre, '_', int2str(session_idx), '_0.mat'];
        post_file = ['../GLM_data/', dataset_post, '/raster_', dataset_post, '_', int2str(session_idx), '_0.mat'];
        raster_pre = load(pre_file, 'rasters').rasters;
        raster_post = load(post_file, "rasters").rasters;
        raster_pre = [raster_pre{:}];
        raster_post = [raster_post{:}];
        spike_pre = load(pre_file, 'spikes').spikes;
        spike_post = load(post_file, "spikes").spikes;
        isi_pre = cellfun(@diff, spike_pre, "UniformOutput", false);
        isi_post = cellfun(@diff, spike_post, "UniformOutput", false);

        N = size(spike_pre, 1);
        %% calc fanos
        % fano = var(isi)/(mean(isi)^2)
        fano_pre = zeros(N, 1);
        fano_post = zeros(N, 1);
        for i=1:N
            isi_pre_i = [isi_pre{i, :}];
            isi_post_i = [isi_post{i, :}];
            fano_pre(i) = var(isi_pre_i)/(mean(isi_pre_i)^2);
            fano_post(i) = var(isi_post_i)/(mean(isi_post_i)^2);
        end

        fano_pre = log(fano_pre);
        fano_post = log(fano_post);
        
        %% plotting
        plot_i = (control_idx-1)*10+session_idx;

        pos1 = plot_i * 3 - 2;
        pos2 = plot_i * 3 - 1;
        group1 = fano_pre;
        group2 = fano_post;
        
        
        for j = 1:length(group1)
            plot([pos1, pos2], [group1(j), group2(j)], 'Color', [0.7 0.7 0.7], 'LineStyle', '-','Marker','none');
        end

        % Plot data points
        scatter(repmat(pos1, size(group1)), group1, 'b', '.', 'MarkerFaceAlpha', 0.3);
        scatter(repmat(pos2, size(group2)), group2, 'r', '.', 'MarkerFaceAlpha', 0.3);

        % Plot boxplots
        boxplot([group1, group2], 'Positions', [pos1, pos2], 'Widths', 0.6);

        % Perform paired t-test (Should change this to ANOVA!)
        [~, p_value] = ttest(group1, group2);
        text(mean([pos1 pos2]), max(max(group1), max(group2)) + 0.5, sprintf('p = %.3f', p_value), 'HorizontalAlignment', 'center');

    end
end


% Customize plot
xticks(1.5:3:43.5);
xticklabels({'Mus 1', 'Mus 2', 'Mus 3', 'Mus 4', 'Mus 5', 'Mus 6', 'Mus 7', 'Mus 8', 'Mus 9', 'Mus 10',...
    'Sal 1', 'Sal 2', 'Sal 3', 'Sal 4', 'Sal 5'});
ylabel('ln(Fano Factor)');
title('Fano Factor');

% Plot dummy points for legend
h1 = scatter(nan, nan, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
h2 = scatter(nan, nan, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
legend([h1, h2], {'Pre', 'Post'}, 'Location', 'best');

hold off;
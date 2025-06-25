models = {'Ctr', 'Sync', 'Mus'};
n_models = length(models);

sig_threshold = 1;

epoch = 1500;
reg = 'L2=2';
kernel = 'DeltaPure';

% f = figure('Visible', 'off', 'Position', [100, 100, 800, 800]);
f = figure('Position', [100, 100, 800, 800]);
t = tiledlayout(3, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

cmap = brewermap(256,'*RdBu');

J_counts = zeros(3, 3, 3, 2); % session, kernel, within/across/all, pos/neg

for model_idx = 1:n_models
    model_name = models{model_idx};
    % model_name_full = ['Simulated_', model_name, '_partial'];
    model_name_full = ['Simulated_', model_name, '_partial_mini'];
    
    % load model
    model_path = ['../GLM_model/', model_name_full, '/GLM_',...
        model_name_full, '_1_', kernel, '_0_', reg, '_', num2str(epoch), '.mat'];
    load(model_path, 'N', 'model_par', 'model_err');

    for k = 1:3
        nexttile;
        J = model_par(:, 1+(k-1)*N+1:1+k*N);
        err = model_err.total(:, 1+(k-1)*N+1:1+k*N);
        imagesc(J);
        clim([-2, 2]);
        colormap(cmap);
        title([model_name, ', Kernel ', num2str(k)]);

        % count significant connections
        % J_within1 = J(1:30, 1:30); J_within2 = J(31:60, 31:60);
        % J_across1 = J(1:30, 31:60); J_across2 = J(31:60, 1:30);
        % err_within1 = err(1:30, 1:30); err_within2 = err(31:60, 31:60);
        % err_across1 = err(1:30, 31:60); err_across2 = err(31:60, 1:30);
        J_within1 = J(1:15, 1:15); J_within2 = J(16:30, 16:30);
        J_across1 = J(1:15, 16:30); J_across2 = J(16:30, 1:15);
        err_within1 = err(1:15, 1:15); err_within2 = err(16:30, 16:30);
        err_across1 = err(1:15, 16:30); err_across2 = err(16:30, 1:15);
        J_counts(model_idx, k, 1, 1) = sum(J_within1(:) > err_within1(:)*sig_threshold)+...
            sum(J_within2(:) > err_within2(:)*sig_threshold);
        J_counts(model_idx, k, 1, 2) = sum(J_within1(:) < -err_within1(:)*sig_threshold)+...
            sum(J_within2(:) < -err_within2(:)*sig_threshold);
        J_counts(model_idx, k, 2, 1) = sum(J_across1(:) > err_across1(:)*sig_threshold)+...
            sum(J_across2(:) > err_across2(:)*sig_threshold);
        J_counts(model_idx, k, 2, 2) = sum(J_across1(:) < -err_across1(:)*sig_threshold)+...
            sum(J_across2(:) < -err_across2(:)*sig_threshold);
        J_counts(model_idx, k, 3, 1) = J_counts(model_idx, k, 1, 1)+J_counts(model_idx, k, 2, 1);
        J_counts(model_idx, k, 3, 2) = J_counts(model_idx, k, 1, 2)+J_counts(model_idx, k, 2, 2);
    end
end

%% bar plot of J counts
f = figure('Position', [100, 100, 800, 800]);
t = tiledlayout(3, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

for model_idx = 1:n_models
    model_name = models{model_idx};
    for k = 1:3
        nexttile;
        bar(squeeze(J_counts(model_idx, k, 1:2, :)), 'grouped', 'BarWidth', 0.8);
        title([model_name, ', Kernel ', num2str(k)]);
        xticks([1, 2]);
        xticklabels({'Within', 'Across', 'All'});
        ylabel('Count');
        legend({'Pos', 'Neg'});
        ylim([0, 200]);
        % text for bars
        for i = 1:2
            for j = 1:2
                text(i+(j-1.5)*0.3, J_counts(model_idx, k, i, j)+20, num2str(J_counts(model_idx, k, i, j)),...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
            end
        end
    end
end
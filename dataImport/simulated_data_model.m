% Generate simulated data for the GLM model
% Validate possible model structures=

% UseGPU = gpuDeviceCount > 0;
UseGPU = false;
partial = true;

for session_idx = 1:10
    %% Generate basic connectivity matrix
    fprintf('Session %d/%d\n', session_idx, 10);
    seed = 137+session_idx; 
    rng(seed);
    
    groups = {'Thal', 'ACC-E', 'ACC-I', 'VLPFC-E', 'VLPFC-I'};
    n_groups = length(groups);
    n_cells = [20, 30, 10, 30, 10];
    n_cells_cumulative = cumsum(n_cells);
    n_cells_cumulative = [0, n_cells_cumulative];
    
    subgroup_size = [10, 10, 10];
    subgroup_size_cumulative = cumsum(subgroup_size);
    subgroup_size_cumulative = [0, subgroup_size_cumulative];
    
    global_J_multi = 1;
    kernel_J_multi = [1, 0.5, 0.25];
    
    % strength parameters
    within_E = 2.7; 
    across_E = 1;
    TC_E = 6; 
    TC_I = 3;
    within_EI = 7.5; 
    within_IE = -5;
    across_EI = 0.5; 
    across_IE = 0;
    
    % baseline h parameters
    thalamus_h = -3.5; % thalamus h
    E_h = -4; % excitatory h
    I_h = -4; % inhibitory h
    
    
    % connection_density = [  0,      0,      0,      0,      0   ;
    %                         0.3,    0.2,    0.15,    0.05,   0   ;
    %                         0.6,    0.2,    0,      0.05,   0   ;
    %                         0.3,    0.05,   0,      0.2,    0.15 ;
    %                         0.6,    0.05,   0,      0.2,    0   ];
    connection_density = [  0,      0,      0,      0,      0   ;
                            0.25,    0.26,    0.4,    0.05,   0   ;
                            0.4,    0.4,    0,      0.05,   0   ;
                            0.25,    0.05,   0,      0.26,    0.4 ;
                            0.4,    0.05,   0,      0.4,    0   ];
    % connection_density = [  0,      0,      0,      0,      0   ;
    %                         0,      0.2,    0.2,    0.05,   0   ;
    %                         0,      0.2,    0,      0.05,   0   ;
    %                         0,      0.05,   0,      0.2,    0.2 ;
    %                         0,      0.05,   0,      0.2,    0   ];
    connection_strength = [ 0,      0,      0,      0,      0   ;
                            TC_E,   within_E,   within_IE,  across_E,   0        ;
                            TC_I,   within_EI,  0,          across_EI,  0        ;
                            TC_E,   across_E,   0,          within_E,   within_IE;
                            TC_I,   across_EI,  0,          within_EI,  0        ];
    
    baseline_h = [thalamus_h, E_h, I_h, E_h, I_h];
    
    N = sum(n_cells);
    J = zeros(N, N, 3);
    h = zeros(N, 1);
    
    
    % generate random J
    for i = 1:N 
        area_i = find(n_cells_cumulative < i, 1, 'last');
    
        if area_i == 2 || area_i == 4
            subgroup_i = find(subgroup_size_cumulative < i-n_cells_cumulative(area_i), 1, 'last');
        else
            subgroup_i = 0;
        end
    
        for j = 1:N
            area_j = find(n_cells_cumulative < j, 1, 'last');
    
            if area_j == 2 || area_j == 4
                subgroup_j = find(subgroup_size_cumulative < j-n_cells_cumulative(area_j), 1, 'last');
            else
                subgroup_j = 0;
            end
    
            for k = 1:3
                if subgroup_i == subgroup_j && area_i == area_j && subgroup_i ~= 0
                    % within subgroup connection
                    density = connection_density(area_i, area_j) * 2.4;
                    strength = connection_strength(area_i, area_j) * 1.6;
                else
                    % between subgroup connection
                    density = connection_density(area_i, area_j);
                    strength = connection_strength(area_i, area_j);
                end
    
                if rand() < density
                    J(i, j, k) = strength*(0.5 + rand());
                end
                % % adjecent neurons have stronger connection
                % if abs(i-j) <= 5 && rand() < 0.3
                %     J(i, j, k) = J(i, j, k) + 0.5 + rand();
                % end
            end
        end
    end
    
    J = J*global_J_multi;
    for k = 1:3
        J(:, :, k) = J(:, :, k) - diag(diag(J(:, :, k))); % remove self connection
        J(:, :, k) = J(:, :, k) * kernel_J_multi(k); % scale by kernel_J_multi
    end
    
    % generate random h
    for i = 1:N
        area_i = find(n_cells_cumulative < i, 1, 'last');
        h(i) = baseline_h(area_i) + randn() - 0.5;
    end
    
    n_trial = 50;
    B = 10000;
    
    % load kernels
    load('../GLM_data/kernel_DeltaPure.mat', 'conn_kernels', 'kernel_len');
    
    conn_kernels_reversed = cell(1, 3);
    for k = 1:3
        conn_kernels_reversed{k} = conn_kernels{k}(:, end:-1:1);
    end
    
    
    %% generate random raster for control, muscimol and sync session
    
    if partial
        session_names = {'Simulated_Ctr_partial', 'Simulated_Sync_partial', 'Simulated_Mus_partial'};
    else
        session_names = {'Simulated_Ctr', 'Simulated_Sync', 'Simulated_Mus'};
    end
    J_ctr = J;
    all_firing_rates = cell(1, 3);
    
    Sync_T = 100; % sync wave period
    thalamus_mask = zeros(1, N);
    thalamus_mask(1:20) = 1;
    all_rasters = cell(1, 3);
    
    for session_type_idx = 1:3
        N = sum(n_cells);
        session_name_full = session_names{session_type_idx};
    
        if strcmp(session_name_full(1:13), 'Simulated_Mus')
            % remove thalamus connections
            J(:, 1:20, :) = 0;
            J_Mus = J;
        elseif strcmp(session_name_full(1:13), 'Simulated_Syn')
            % % uniform thalamocortical connections
            % sync_J = mean(J(21:100, 1:20, :), 2);
            % J(21:100, 1:20, :) = repmat(sync_J, 1, 20, 1);
            J_Sync = J;
        end
    
        rasters = cell(1, n_trial);
        trial_len = ones(1, n_trial)*B;
        spikes = cell(1, n_trial);
        firing_rates = zeros(N, n_trial);
    
        % GPU if available
        if UseGPU
            fprintf("Using GPU\n");
            J = gpuArray(J);
            h = gpuArray(h);
            conn_kernels_reversed = cellfun(@gpuArray, conn_kernels_reversed, 'UniformOutput', false);
        end
    
        parfor trial_idx = 1:n_trial
            fprintf("%s, trial %d\n", session_name_full, trial_idx);
            raster = zeros(N, B);
    
            if UseGPU
                raster = gpuArray(raster);
            end
    
            for t = 1:B
                htot = h;
    
                if strcmp(session_name_full(1:13), 'Simulated_Syn')
                    % sync wave
                    if mod(t, Sync_T) < Sync_T/8
                        htot = htot + 1.5*thalamus_mask';
                    else
                        htot = htot - 0.6*thalamus_mask';
                    end
                end
    
                if t>=kernel_len
                    for k = 1:3
                        kernel = conn_kernels_reversed{k};
                        % raster_conv = raster(:, t-kernel_len+1:t);
                        % convolved = sum(raster_conv .* kernel, 2);
                        convolved = sum(raster(:, t-kernel_len+1:t) .* kernel, 2);
                        h_k = J(:, :, k)*convolved;
                        htot = htot + h_k;
                    end
                end
                lambda = exp(htot);
                p = lambda./(1+lambda);
                raster(:, t) = rand(N, 1) < p;
            end
    
            if UseGPU
                raster = gather(raster);
            end
    
            rasters{trial_idx} = raster;
            firing_rates(:, trial_idx) = sum(raster, 2)/B;
        end
    
        all_rasters{session_type_idx} = rasters;
        all_firing_rates{session_type_idx} = firing_rates;
    
        if partial
            % remove thalamus and inhibitory neurons
            for i=1:n_trial
                rasters{i}=[rasters{i}(21:50, :); rasters{i}(61:90, :)];
            end
            firing_rates = [firing_rates(21:50, :); firing_rates(61:90, :)];
            N = 60;
        end
    
        cell_id = cell(1, N);
        cuetype = zeros(1, n_trial);
        channel = 1:N;
        % borders = [-0.5];
    
        check_path(['../GLM_data/', session_name_full]);
        save_file = ['../GLM_data/', session_name_full, '/raster_', session_name_full, '_', num2str(session_idx), '_0.mat'];
        save(save_file, 'N', 'rasters', 'J', 'h', "n_trial", "session_name_full", "trial_len", "channel", "firing_rates");
    
        % borders = [20.5, 60.5];
        borders = [30.5, 30.5];
        save(['../GLM_data/', session_name_full, '/borders_', session_name_full, '_', num2str(session_idx), '.mat'], "borders");
    
    end
    
    N = sum(n_cells);
    
    model_names = {'Ctr', 'Sync', 'Mus'};
    % plot J matrix
    cmap = brewermap(256,'*RdBu');
    all_J = {J_ctr, J_Sync, J_Mus};
    figure('Position', [100, 100, 800, 800]);
    t = tiledlayout(3, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    for i = 1:3
        for k = 1:3
            model_name = model_names{i};
            nexttile;
            hold on;
            imagesc(squeeze(all_J{i}(:, :, k)));
            % title(['J', num2str(k)]);
            colorbar;
            colormap(cmap);
            clim([-7.5, 7.5])
            axis square;
            set(gca, 'XTick', [], 'YTick', []);
            title([model_name, ', Kernel ', num2str(k)]);

            % area borders
            for j = 1:(n_groups+1)
                plot([n_cells_cumulative(j), n_cells_cumulative(j)]+0.5, [0, N]+0.5, 'k', 'LineWidth', 1);
                plot([0, N]+0.5, [n_cells_cumulative(j), n_cells_cumulative(j)]+0.5, 'k', 'LineWidth', 1);
            end

            % reverse y axis
            set(gca, 'YDir', 'reverse');
            xlim([0.5, N+0.5]);
            ylim([0.5, N+0.5]);

            hold off;
        end
    end

    % plot connection counts
    lim = [0.5, 100.5];
    borders = [0.5, 20.5, 50.5, 60.5, 90.5, 100.5];
    areas = [20, 50; 60, 90];
    all_J = {J_ctr, J_Sync, J_Mus};
    figure('Position', [100, 100, 800, 800]);
    t = tiledlayout(3, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    for session_i = 1:3
        for k = 1:3
            nexttile;
            hold on;
            J = squeeze(all_J{session_i}(:, :, k));

            % count connections
            J_within1 = J((areas(1, 1)+1):areas(1, 2), (areas(1, 1)+1):areas(1, 2));
            J_within2 = J((areas(2, 1)+1):areas(2, 2), (areas(2, 1)+1):areas(2, 2));
            J_across1 = J((areas(1, 1)+1):areas(1, 2), (areas(2, 1)+1):areas(2, 2));
            J_across2 = J((areas(2, 1)+1):areas(2, 2), (areas(1, 1)+1):areas(1, 2));
            err = 1e-5;

            within_pos = nnz(J_within1(:) > err)+nnz(J_within2(:) > err);
            within_neg = nnz(J_within1(:) < -err)+nnz(J_within2(:) < -err);
            across_pos = nnz(J_across1(:) > err)+nnz(J_across2(:) > err);
            across_neg = nnz(J_across1(:) < -err)+nnz(J_across2(:) < -err);

            J_counts = [within_pos, within_neg; across_pos, across_neg; within_pos+across_pos, within_neg+across_neg];

            b = bar([1, 2], J_counts(1:2, :), 'grouped', 'BarWidth', 0.8,'FaceColor', 'flat');
            % set colors for bars: red for positive, blue for negative
            b(1).CData(1, :) = [0.9, 0, 0]; % red for positive
            b(2).CData(1, :) = [0, 0, 0.9]; % blue for negative
            b(1).CData(2, :) = [0.9, 0, 0]; % red for positive
            b(2).CData(2, :) = [0, 0, 0.9]; % blue for negative

            title([model_name, ', Kernel ', num2str(k)]);
            xticks([1, 2]);
            xticklabels({'Within', 'Across', 'All'});
            ylabel('Count');
            legend({'Pos', 'Neg'});
            ylim([0, 1200]);
            % text for bars
            for i = 1:2
                for j = 1:2
                    text(i+(j-1.5)*0.3, J_counts(i, j)+20, num2str(J_counts(i, j)),...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
                end
            end

            hold off;
        end
    end

    % Plot firing rates
    figure('Position', [100, 100, 800, 800]);
    hold on;
    for i = 1:3
        plot(mean(all_firing_rates{i}, 2), 'o-');
    end
    hold off;
    title('Firing rates');
    xlabel('Neuron ID');
    ylabel('Firing rate (Hz)');
    legend(session_names, 'Location', 'Best');

    % Plot raster
    figure('Position', [100, 100, 800, 800]);
    t = tiledlayout(3, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    for i = 1:3
        rasters = all_rasters{i};
        nexttile;
        % plot thalamus raster in first 10 trials, first 10 neurons
        concatenated_raster = zeros(10, B);
        shift = 0:2:18;
        for j = 1:1
            concatenated_raster(:, (j-1)*B+1:j*B) = rasters{j}(1:10, :);
        end
        imagesc(concatenated_raster);
        title(['Session ', num2str(i)]);
        xlim();
    end 

end

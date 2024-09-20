%% load data
session_num = 10;
kernel_name='expDecay10';
reg.name='L1=5';
epoch=2500;
shuffle_range=1:4;
shuffle_size = length(shuffle_range);

J_mean = ones(7, session_num, 2, 3) * NaN; % conn_type, session, sub_idx, neg/pos/abs
J_err = ones(7, session_num, 2, 3) * NaN;
sig = zeros(7, session_num, 3); % conn_type, session, neg/pos/abs

for session = 1:session_num
    %% load models and align pre-post
    taskname = 'eyeOpen';
    modelfile_pre = ['../GLM_model/', taskname, '/GLM_', taskname, '_', int2str(session), ...
        '_', kernel_name, '_0_', reg.name, '_', int2str(epoch), '.mat'];
    taskname = 'eyeClose';
    modelfile_post = ['../GLM_model/', taskname, '/GLM_', taskname, '_', int2str(session), ...
        '_', kernel_name, '_0_', reg.name, '_', int2str(epoch), '.mat'];

    load(modelfile_pre, "model_par", "N");
    model_par_pre = model_par;
    N_pre=N;
    load(modelfile_post, "model_par", "N");
    model_par_post = model_par;
    N_post=N;

    taskname = 'eyeOpen';
    borderfile = ['../GLM_data/', taskname,'/borders_', taskname, '_', ...
    int2str(session),'.mat'];
    load(borderfile, "borders");
    sortfile = ['../GLM_data/', taskname,'/sortidx_', taskname, '_', ...
    int2str(session), '_', kernel_name, '.mat'];
    load(sortfile, "sort_idx");

    % round borders to integers
    A_T_border = borders(1) + 0.5;
    T_P_border = borders(2) + 0.5;
    borders = [1, round(borders+0.5, 0)];

    J_pre = model_par_pre(:, 2:N_pre+1);
    J_post = model_par_post(:, 2:N_post+1);

    % load shuffled data
    J_sfl_pre = zeros([size(J_pre), shuffle_size]);
    J_sfl_post = zeros([size(J_post), shuffle_size]);

    for i=1:shuffle_size
        taskname = 'eyeOpen';
        path_sfl_pre = ['../GLM_model/', taskname, '/GLM_', taskname, '_', int2str(session), ...
        '_', kernel_name, '_', int2str(i), '_', reg.name, '_', int2str(epoch), '.mat'];
        par_sfl = load(path_sfl_pre).model_par;
        J_sfl_pre(:, :, i) = par_sfl(:, 2:N_pre+1);

        taskname = 'MuscimolPost_cortex';
        path_sfl_post = ['../GLM_model/', taskname, '/GLM_', taskname, '_', int2str(session), ...
        '_', kernel_name, '_', int2str(i), '_', reg.name, '_', int2str(epoch), '.mat'];
        par_sfl = load(path_sfl_post).model_par;
        J_sfl_post(:, :, i) = par_sfl(:, 2:N_post+1);
    end

    % statistical test
    significant_pre = zeros(size(J_pre));
    significant_post = zeros(size(J_post));
    
    for i=1:N_pre
        for j=1:N_pre
            % one-sample t-test
            h = ttest(reshape(J_sfl_pre(i, j, :)-J_pre(i, j), [1, shuffle_size]), ...
                0, "Alpha", 0.001);
            significant_pre(i, j) = h;
    
        end
    end

    for i=1:N_post
        for j=1:N_post
            % one-sample t-test
            h = ttest(reshape(J_sfl_post(i, j, :)-J_post(i, j), [1, shuffle_size]), ...
                0, "Alpha", 0.001);
            significant_post(i, j) = h;
    
        end
    end

    J_pre = (J_pre - mean(J_sfl_pre, 3)).*significant_pre;
    J_post = (J_post - mean(J_sfl_post, 3)).*significant_post;

    % fill NaN to post thalamus area
    assert(N_pre-N_post == borders(3)-borders(2));
    J_post = [J_post(:, 1:borders(2)-1), ones(N_post, N_pre-N_post)*NaN, J_post(:, borders(2):end)];
    J_post = [J_post(1:borders(2)-1, :); ones(N_pre-N_post, N_pre)*NaN; J_post(borders(2):end, :)];
    N = N_pre;

    % sorting
    J_pre = J_pre(sort_idx, sort_idx);
    J_post = J_post(sort_idx, sort_idx);
   
    %% assign groups
    Js = {J_pre, J_post};
    conn_groups = cell(1, 2);
    mask_abs = ((J_pre~=0)|(J_post~=0)) & ~isnan(J_pre) & ~isnan(J_post);
    for sub_idx=1:2
        J = Js{sub_idx};
        conn_group = zeros(N);
        mask_neg = (-J > 0) & ~isnan(J_pre) & ~isnan(J_post);
        mask_pos = (J > 0) & ~isnan(J_pre) & ~isnan(J_post);
    
        for i = 1:N
            for j = 1:N
                if i < A_T_border && j < A_T_border
                    conn_group(i, j) = 1; % ACC-ACC, R
                end
                if i >= A_T_border && i < T_P_border && j >= A_T_border && j < T_P_border
                    conn_group(i, j) = 2; % Thalamus-Thalamus, G
                end
                if i >= T_P_border && j >= T_P_border
                    conn_group(i, j) = 3; % VLPFC-VLPFC, B
                end
                if i < A_T_border && j >= A_T_border && j < T_P_border
                    conn_group(i, j) = 4; % ACC-Thalamus, Y
                    conn_group(j, i) = 4; 
                end
                if i < A_T_border && j >= T_P_border
                    conn_group(i, j) = 5; % ACC-VLPFC, M
                    conn_group(j, i) = 5; 
                end
                if i >= A_T_border && i < T_P_border  && j >= T_P_border
                    conn_group(i, j) = 6; % Thalamus-VLPFC, C
                    conn_group(j, i) = 6; 
                end
            end
        end
        conn_groups{sub_idx}=conn_group;

        % calculate mean J and err
        for i = 1:6
            J_mean(i, session, sub_idx, 1) = mean(J((conn_group==i) & ~isnan(J) & mask_neg));
            J_mean(i, session, sub_idx, 2) = mean(J((conn_group==i) & ~isnan(J) & mask_pos));
            J_mean(i, session, sub_idx, 3) = mean(abs(J((conn_group==i) & ~isnan(J) & mask_abs)));
            J_err(i, session, sub_idx, 1) = std(J((conn_group==i) & ~isnan(J) & mask_neg))/sqrt(length(J((conn_group==i) & ~isnan(J) & mask_neg)));
            J_err(i, session, sub_idx, 2) = std(J((conn_group==i) & ~isnan(J) & mask_pos))/sqrt(length(J((conn_group==i) & ~isnan(J) & mask_pos)));
            J_err(i, session, sub_idx, 3) = std(abs(J((conn_group==i) &  ~isnan(J) & mask_abs)))/sqrt(length(abs(J((conn_group==i) &  ~isnan(J) & mask_abs))));
        end

        J_mean(7, session, sub_idx, 1) = mean(J((conn_group~=0) & ~isnan(J) & mask_neg));
        J_mean(7, session, sub_idx, 2) = mean(J((conn_group~=0) & ~isnan(J) & mask_pos));
        J_mean(7, session, sub_idx, 3) = mean(abs(J((conn_group~=0) & ~isnan(J) & mask_abs)));
        J_err(7, session, sub_idx, 1) = std(J((conn_group~=0) & ~isnan(J) & mask_neg))/sqrt(length(J((conn_group~=0) & ~isnan(J) & mask_neg)));
        J_err(7, session, sub_idx, 2) = std(J((conn_group~=0) & ~isnan(J) & mask_pos))/sqrt(length(J((conn_group~=0) & ~isnan(J) & mask_pos)));
        J_err(7, session, sub_idx, 3) = std(abs(J((conn_group~=0) & ~isnan(J) & mask_abs)))/sqrt(length(abs(J((conn_group~=0) &  ~isnan(J) & mask_abs))));

        % fprintf("Session%d %s, mean: %.3f\n", ...
        %     session, sub_names{sub_idx}, mean(J(~isnan(J))));
    end

end

%% Plot 1 pos and neg

figure('Visible', 'off', 'Position', [100, 100, 1600, 1200]);
palette = ["red", "blue"];
titles = ["ACC-ACC", "Tha-Tha", "PFC-PFC", "ACC-Tha", "ACC-PFC", "Tha-PFC", "All"];
plot_pos = [1, 5, 9, 2, 3, 6, 7];
x=1:(session_num+1);
for i=1:7
    subplot(3, 3, plot_pos(i));
    hold on
    for sub_idx = 1:2
    % set(gcf, 'color', [alpha_val alpha_val alpha_val]);
        y = J_mean(i, :, sub_idx, 1);
        y(session_num+1) = mean(y, "omitnan");
        y_err = J_err(i, :, sub_idx, 1);
        y_err(session_num+1) = std(y(1:session_num), "omitnan");
        style = '-';
        mark = 'o';
        if sub_idx == 2
            style = '--';
            mark = 'x';
        end
        errorbar(x(~isnan(y))+0.2*(sub_idx-1.5), y(~isnan(y)), y_err(~isnan(y)),'Marker', mark, 'LineStyle', style, "Color", palette(sub_idx));
        % plot(session_num+1, mean(y(~isnan(y))), 'Marker', mark, "Color", palette(sub_idx));

        y = J_mean(i, :, sub_idx, 2);
        y(session_num+1) = mean(y, "omitnan");
        y_err = J_err(i, :, sub_idx, 1);
        y_err(session_num+1) = std(y(1:session_num), "omitnan");
        style = '-';
        mark = 'o';
        if sub_idx == 2
            style = '--';
            mark = 'x';
        end
        errorbar(x(~isnan(y))+0.2*(sub_idx-1.5), y(~isnan(y)), y_err(~isnan(y)),'Marker', mark, 'LineStyle', style, "Color", palette(sub_idx));
        % plot(session_num+1, mean(y(~isnan(y))), 'Marker', mark, "Color", palette(sub_idx));
    end
    xlabel("session number");
    ylabel("mean J");
    title(titles(i));
    if i==7
        legend(["pre", "", "post", ""], "FontSize", 8);
    end
    xlim([0, session_num+2]);
    xticks(1:(session_num+1));
    xticklabels(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "mean"]);
    set(gca, 'XAxisLocation', 'origin');
    hold off
end
saveas(gcf, sprintf('../figures/J_all_posinega.png'));
close gcf;

%% Plot 2 histogram of abs J

figure('Visible', 'off', 'Position', [100, 100, 1600, 1200]);
palette = ["red", "blue"];
titles = ["ACC-ACC", "Tha-Tha", "PFC-PFC", "ACC-Tha", "ACC-PFC", "Tha-PFC", "All"];
plot_pos = [1, 5, 9, 2, 3, 6, 7];
x=1:(session_num+1);
for i=1:7
    subplot(3, 3, plot_pos(i));
    hold on
    for sub_idx = 1:2
    % set(gcf, 'color', [alpha_val alpha_val alpha_val]);
        y = J_mean(i, :, sub_idx, 3);
        y(session_num+1) = mean(y, "omitnan");
        y_err = J_err(i, :, sub_idx, 3);
        y_err(session_num+1) = std(y(1:session_num), "omitnan");
        bar(x(~isnan(y))+0.2*(sub_idx-1.5), y(~isnan(y)), 0.2, "FaceColor", palette(sub_idx));
        errorbar(x(~isnan(y))+0.2*(sub_idx-1.5), y(~isnan(y)), y_err(~isnan(y)), 'LineStyle', 'none', "Color", [0 0 0]);
    end
    xlabel("session number");
    ylabel("mean J");
    title(titles(i));
    if i==7
        legend(["pre", "", "post", ""], "FontSize", 8);
    end
    xlim([0, session_num+2]);
    xticks(1:(session_num+1));
    xticklabels(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "mean"]);
    % set(gca, 'XAxisLocation', 'origin');
    hold off
end
saveas(gcf, sprintf('../figures/J_all_abs.png'));
close gcf;

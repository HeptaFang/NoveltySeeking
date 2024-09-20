% inward curve and outward curve, for 3 areas separetely

% load data
kernel_name='expDecay10';
reg.name='L1=5';
epoch=2500;

for session = 1:10
    fig = figure("Visible", "off");
    set(fig, 'PaperPosition', [0, 0, 20, 6]);
    tiles = tiledlayout(2, 3);

    taskname = 'MuscimolPre_full';
    modelfile_pre = ['../GLM_model/', taskname, '/GLM_', taskname, '_', int2str(session), ...
        '_', kernel_name, '_0_', reg.name, '_', int2str(epoch), '.mat'];
    taskname = 'MuscimolPost_cortex';
    modelfile_post = ['../GLM_model/', taskname, '/GLM_', taskname, '_', int2str(session), ...
        '_', kernel_name, '_0_', reg.name, '_', int2str(epoch), '.mat'];

    load(modelfile_pre, "model_par", "N");
    model_par_pre = model_par;
    N_pre=N;
    load(modelfile_post, "model_par", "N");
    model_par_post = model_par;
    N_post=N;

    taskname = 'MuscimolPre_full';
    borderfile = ['../GLM_data/', taskname,'/borders_', taskname, '_', ...
    int2str(session),'.mat'];
    load(borderfile, "borders");
    sortfile = ['../GLM_data/', taskname,'/sortidx_', taskname, '_', ...
    int2str(session), '_', kernel_name, '.mat'];
    load(sortfile, "sort_idx");

    % round borders to integers
    borders = [1, round(borders+0.5, 0)];

    J_pre = model_par_pre(:, 2:N_pre+1);
    J_post = model_par_post(:, 2:N_post+1);

    % fill NaN to post thalamus area
    assert(N_pre-N_post == borders(3)-borders(2));
    J_post = [J_post(:, 1:borders(2)-1), ones(N_post, N_pre-N_post)*NaN, J_post(:, borders(2):end)];
    J_post = [J_post(1:borders(2)-1, :); ones(N_pre-N_post, N_pre)*NaN; J_post(borders(2):end, :)];

    % sorting
    J_pre = J_pre(sort_idx, sort_idx);
    J_post = J_post(sort_idx, sort_idx);

    plot_colors = {[1 0 0.25], [1 0.75 0], [0, 0, 1]};
    plot_names = {'ACC', 'Thalamus', 'VLPFC'};
    
    for i=1:3
        area_s = borders(i);
        area_e = borders(i+1)-1;
        inward_pre = mean(J_pre(:, area_s:area_e), 2, "omitnan");
        outward_pre = mean(J_pre(area_s:area_e, :), 1, "omitnan");
        inward_post = mean(J_post(:, area_s:area_e), 2, "omitnan");
        outward_post = mean(J_post(area_s:area_e, :), 1, "omitnan");

        % inward plot
        nexttile(i);
        hold on
        % plot(inward_pre, "LineStyle", "--", "Color",plot_colors{i}, "Marker",".", "DisplayName",[plot_names{i}, ' pre']);
        % plot(inward_post, "LineStyle", "-", "Color",plot_colors{i}, "Marker",".", "DisplayName",[plot_names{i}, ' post']);
        plot(inward_pre, "b.-", "DisplayName", 'pre');
        plot(inward_post, "r.-", "DisplayName", 'post');
        title(['cell i to ', plot_names{i}]);
        legend("Location","bestoutside");
        xlabel("From cell idx");
        ylabel("Ave inward J");
        set(gca, "XAxisLocation", "origin");
        xline(borders(2:3)-0.5, "HandleVisibility","off");
        hold off

        % outward plot
        nexttile(3+i);
        hold on
        % plot(outward_pre, "LineStyle", "--", "Color",plot_colors{i}, "Marker",".", "DisplayName",[plot_names{i}, ' pre']);
        % plot(outward_post, "LineStyle", "-", "Color",plot_colors{i}, "Marker",".", "DisplayName",[plot_names{i}, ' post']);
        plot(outward_pre, "b.-", "DisplayName",'pre');
        plot(outward_post, "r.-", "DisplayName",'post');
        title([plot_names{i}, ' to cell i']);
        legend("Location","bestoutside");
        xlabel("To cell idx");
        ylabel("Ave outward J");
        set(gca, "XAxisLocation", "origin");
        xline(borders(2:3)-0.5, "HandleVisibility","off");
        hold off
    end

    % suptitle(['Session ', int2str(session)]);

    fig_path = ['../figures/GLM/', 'Jmatrix'];
    check_path(fig_path);
    fig_file = [fig_path, '/Jmat_session', int2str(session), '_',...
            kernel_name, '_', reg.name, '_', int2str(epoch), '.png'];
    % exportgraphics(fig,fig_file);
    print(fig, fig_file,'-dpng', '-r300');

end
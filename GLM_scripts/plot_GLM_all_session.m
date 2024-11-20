
dataset_name = 'MuscimolPre_full'
session_num = 10;
kernel_name = 'expDecay10';
reg.name = 'L1=5';
epoch = 2500;

% figure
fig = figure("Visible","off","Position",[0, 0, 3000, 3000]);
% fig = figure;


for session = 1:session_num
    % load data
    data_path = ['../GLM_data/', dataset_name, '/GLMdata_', dataset_name, '_', ...
        int2str(session),'_', kernel_name, '_0.mat'];
    model_path=['../GLM_model/', dataset_name, '/GLM_', dataset_name, '_', ...
        int2str(session), '_', kernel_name, '_0_', ...
        reg.name, '_', int2str(epoch)];
    border_path=['../GLM_data/', dataset_name,'/borders_', dataset_name, '_', ...
        int2str(session),'.mat'];
    sort_path = ['../GLM_data/', dataset_name, '/sortidx_', dataset_name, '_', ...
        int2str(session),'_', kernel_name, '.mat'];

    load(data_path, 'raster');
    load(model_path, "model_par", "N");
    load(border_path, "borders");
    load(sort_path, "sort_idx");

    % matrix and firing rate
    firing_rate = mean(raster, 2)*1000;
    J = model_par(:, 2:N+1);

    firing_rate = firing_rate(sort_idx);
    J = J(sort_idx, sort_idx);

    % plot
    if session<=5
        pos=session;
    else
        pos=session+5;
    end

    subplot(4,5, pos);
    plot(1:N, firing_rate, '.-');
    title(['session', int2str(session)]);
    xlabel('cell idx');
    ylabel('firing rate (Hz)')
    for b=borders
        xline(b);
    end
    xlim([0.5, N+0.5]);
    ylim([0, 40]);
    axis square;

    subplot(4,5, pos+5);
    imagesc(J);
    clim([-1, 1]);
    cmap=brewermap(256,'*RdBu');
    colormap(gca, cmap);
    % colorbar;
    title('J');
    xlabel('from idx');
    ylabel('to idx');
    for b=borders
        line([0.5, N+0.5], [b, b], 'Color', 'k');
        line([b, b], [0.5, N+0.5], 'Color', 'k');
    end
    axis square;

end

sgtitle(dataset_name);

fig_path = ['../figures/GLM/', dataset_name, '/connectivity_' dataset_name, '.png'];
saveas(fig, fig_path);
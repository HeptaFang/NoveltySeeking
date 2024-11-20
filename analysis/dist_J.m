% inward curve and outward curve, for 3 areas separetely

% load data
kernel_name='expDecay10';
reg.name='L1=5';
epoch=2500;

% Initialize figure
fig = figure("Visible", "off");
set(fig, 'PaperPosition', [0, 0, 15, 10]);

% Colors for pre and post states
pre_color = 'b';
post_color = 'r';

% Define subplot positions
subplot_positions = [1, 2, 3; 4, 2, 6; 7, 8, 3];

for session = 1:10

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
    channelfile = ['../GLM_data/', taskname, '/raster_', taskname, '_', int2str(session), ...
    '_0.mat'];
    load(channelfile, "channel", "N");

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
    channel = channel(sort_idx);

    % plot_colors = {[1 0 0.25], [1 0.75 0], [0, 0, 1]};
    plot_names = {'ACC', 'Thalamus', 'VLPFC'};
    
    % Iterate through the areas
    for area_i = 1:3
        for area_j = 1:3
            if area_j~=area_i
                continue;
            end
            % Initialize variables to store distances and connectivities
            distances_pre = [];
            connectivities_pre = [];
            distances_post = [];
            connectivities_post = [];
            
            % Get the unit indices for areas
            units_i = borders(area_i):borders(area_i+1)-1;
            units_j = borders(area_j):borders(area_j+1)-1;
            
            % Iterate through the units in the current areas
            for unit_i = units_i
                for unit_j = units_j
                    if unit_i ~= unit_j
                        % Calculate distance
                        distance = abs(channel(unit_i) - channel(unit_j));
                        
                        % Get connectivity strengths
                        connectivity_pre = J_pre(unit_i, unit_j);
                        connectivity_post = J_post(unit_i, unit_j);
                        
                        % Store distances and connectivities
                        distances_pre = [distances_pre, distance];
                        connectivities_pre = [connectivities_pre, connectivity_pre];
                        distances_post = [distances_post, distance];
                        connectivities_post = [connectivities_post, connectivity_post];
                    end
                end
            end
            
            % Plot the results in the appropriate subplot
            subplot(1, 3, subplot_positions(area_i, area_j));
            hold on;
            scatter(distances_pre, connectivities_pre, 5, '.', pre_color);
            scatter(distances_post, connectivities_post, 5, '.', post_color);
            xlabel('Distance');
            ylabel('Connectivity');
            ylim([-2, 2]);
            title(sprintf('%s to %s', plot_names{area_j}, plot_names{area_i}));
            legend('Pre', 'Post');
            hold off;
        end
    end



end

sgtitle('Connectivity vs Distance');
fig_path = ['../figures/GLM/', 'Jmatrix'];
check_path(fig_path);
fig_file = [fig_path, '/Jmat_dist', int2str(session), '_',...
        kernel_name, '_', reg.name, '_', int2str(epoch), '.png'];
% exportgraphics(fig,fig_file);
print(fig, fig_file,'-dpng', '-r300');

% Adjust layout
% saveas(gcf, 'connectivity_vs_distance.png');
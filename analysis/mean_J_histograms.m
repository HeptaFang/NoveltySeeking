function mean_J_histograms(session_type)
root_path = '../';
kernel = 'Delta';
reg = 'L2=3e3';
epoch = '2500';
area_names = {'ACC', 'Thalamus', 'VLPFC'};
% load data
states = {'Decision', 'Info', 'InfoAnti', 'InfoResp', 'RestClose', 'RestOpen'};
n_states = length(states);
if strcmp(session_type, 'Muscimol')
    n_session = 10;
else
    n_session = 5;
end

% % load metadata
% if strcmp(prepost, 'Pre')
%     session_stage_full = [session_type, 'Pre', states{1}, '_full'];
% else
%     session_stage_full = [session_type, 'Post', states{1}, '_cortex'];
% end
% 
% % load model
% file_path = [root_path, 'GLM_model/', session_stage_full,...
%         '/GLM_', session_stage_full, '_', '1', '_',...
%         kernel, '_0_', reg, '_', epoch, '.mat'];
% load(file_path,"n_conn_kernel", "kernel_len", "N");
n_conn_kernel = 3;
J_data = cell(3, 3, n_conn_kernel, n_states, n_session, 2);

all_J = [];
J_state = [];
J_session = [];
J_kernel = [];
J_i = [];
J_j = [];

prepost_all = {'Pre', 'Post'};
for prepost_idx = 1:2
    prepost = prepost_all{prepost_idx};
    for session_idx = 1:n_session
        fprintf('Loading session %d\n', session_idx);
        for state_idx = 1:n_states
            state = states{state_idx};
            if strcmp(prepost, 'Pre')
                session_stage_full = [session_type, 'Pre', state, '_full'];
            else
                session_stage_full = [session_type, 'Post', state, '_cortex'];
            end
    
            file_path = [root_path, 'GLM_model/', session_stage_full,...
                '/GLM_', session_stage_full, '_', num2str(session_idx), '_',...
                kernel, '_0_', reg, '_', epoch, '.mat'];
    
            fprintf('Loading %s\n', session_stage_full);
    
            load(file_path, "model_par", "n_PS_kernel", "kernel_len", "N");
            load([root_path, 'GLM_data/', session_stage_full,'/borders_', session_stage_full, '_', ...
                    num2str(session_idx),'.mat'], "borders");
            
            borders = [1, borders+0.5]; % area i is from borders(i) to borders(i+1)
            n_area = length(borders) - 1;
            % if strcmp(prepost, 'Post')
            %     assert(n_area == 2, 'Only 2 areas in Post sessions');
            % else
            %     assert(n_area == 3, 'Only 3 areas in Pre sessions');
            % end
    
            for kernel_idx = 1:n_conn_kernel
                J_mat = model_par(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));
                for i = 1:n_area
                    for j = 1:n_area
                        if strcmp(prepost, 'Post') && i == 2
                            i_eff = 3;
                        else
                            i_eff = i;
                        end
                        if strcmp(prepost, 'Post') && j == 2
                            j_eff = 3;
                        else
                            j_eff = j;
                        end

                        J_area = J_mat(borders(i):borders(i+1)-1, borders(j):borders(j+1)-1);
                        J_data{i_eff, j_eff, kernel_idx, state_idx, session_idx, prepost_idx} = J_area(:);
                        all_J = [all_J; J_area(:)];
                        J_state = [J_state, repelem(state_idx, numel(J_area))];
                        J_session = [J_session, repelem(session_idx, numel(J_area))];
                        J_kernel = [J_kernel, repelem(kernel_idx, numel(J_area))];
                        J_i = [J_i, repelem(i_eff, numel(J_area))];
                        J_j = [J_j, repelem(j_eff, numel(J_area))];
                    end
                end
            end
        end
    end
end

% if strcmp(prepost, 'Pre')
%     session_stage_full = [session_type, 'Pre', '_full'];
% else
%     session_stage_full = [session_type, 'Post', '_cortex'];
% end
n_area = 3;
ylim_all_kernel = {[-0.01, 0.25], [-0.01, 0.15], [-0.01, 0.12]};
for kernel_idx = 1:n_conn_kernel
    % Each kernel is one figure
    f = figure("Visible", "off","Position",[0, 0, 1600, 900]);
    % ij are subplots
    tiles = tiledlayout(n_area, n_area);

    for i = 1:n_area
        for j = 1:n_area
            fprintf('Kernel%d, i=%d, j=%d\n', kernel_idx, i, j);
            nexttile;
            % % violin plot
            % filter = J_i == i & J_j == j & J_kernel == kernel_idx;
            % data = all_J(filter);
            % data_state = J_state(filter);
            % data_session = J_session(filter);
            % data_by_session = cell(1, n_session);
            % for session_idx = 1:n_session
            %     data_by_session{session_idx} = data(data_session == session_idx);
            %     if size(data_by_session{session_idx}, 1) == 0
            %         data_by_session{session_idx} = 0;
            %     end
            % end
            % violinplot(data_session, data, GroupByColor=data_state);
            % % daviolinplot(data_by_session,'groups',data_state,'outsymbol','k+',...
            % %     'xtlabels', 1:n_session,'scatter',2,'jitter',1,...
            % %     'box',1,'boxcolors','same', 'color', c,'scattercolors','same',...
            % %     'boxspacing',1.1,'legend',states);
            % % legend(states);

            % barplot
            means = zeros(2, n_states);
            errors = zeros(2, n_states);
            for prepost_idx = 1:2
                for state_idx = 1:n_states
                    data = [];
                    for session_idx = 1:n_session
                        data = [data; J_data{i, j, kernel_idx, state_idx, session_idx, prepost_idx}(:)];
                    end
                    means(prepost_idx, state_idx) = mean(data);
                    errors(prepost_idx, state_idx) = std(data)/sqrt(length(data));
                end
            end
            for state_idx = 1:n_states
                shift = (state_idx-3.5)*0.12;
                hold on;
                bar((1:2)+shift, means(:, state_idx), 0.12, '');
                errorbar((1:2)+shift, means(:, state_idx), errors(:, state_idx), 'LineStyle', 'none', "Color", [0 0 0],"CapSize",2);
                hold off;
            end
            title([area_names{j}, ' to ', area_names{i}]);
            xticks(1:2);
            xticklabels({'Pre', 'Post'});
            ylabel('Average J');
            if i==2 && j==2
                % state_legends = {'Decision', '', 'Info', '',...
                %     'InfoAnti', '', 'InfoResp', '',...
                %     'RestClose', '', 'RestOpen', ''};
                state_legends = {'before choice', '', 'choice to reward', '',...
                    'choice to infocue', '', 'infocue to reward', '',...
                    'eye close', '', 'eye open', ''};
                legend(state_legends);
            end
            % ylim(ylim_all_kernel{kernel_idx});
            % ylim(ylim_all_kernel{1});
        end
    end
    title(tiles, [session_type,' Kernel ', num2str(kernel_idx)]);

    fig_folder = [root_path, 'figures/GLM/', session_type];
    check_path(fig_folder);
    saveas(f, [fig_folder, '/J_histograms_kernel_', num2str(kernel_idx), '.png']);
end
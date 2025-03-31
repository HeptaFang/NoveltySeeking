% session_type = 'Muscimol';
session_type = 'Saline';
root_path = '../';
kernel = 'DeltaPure';
reg = 'L2=2';
epoch = '2500';
use_filter = false;
filter_threshold = 1;

% area_names = {'ACC', 'Thalamus', 'VLPFC'};
area_names = {'ACC', 'VLPFC'};
% load data
% states = {'Offer1', 'Offer2', 'Decision', 'InfoAnti', 'InfoResp', 'Reward', 'RandomA', 'RandomB'};
% states = {'Offer1', 'Offer2', 'Decision', 'RandomShort', 'RandomLong', 'RandomA', 'RandomB', 'RestOpen', 'RestClose'};
% states = {'RandomA', 'RandomShort', 'RandomLong'};
% states = {'Simulated', 'Simulated_higher', 'RandomLong'};
states = {'Task', 'RestOpen', 'RestClose'};
state_names = {'Task', 'Eyes\_Open', 'Eyes\_Closed'};
aligns = {'AlignFirst', 'AlignLast', 'AlignRandom'};
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
J_data = cell(3, 3, n_conn_kernel, n_states, 3, n_session, 3); % (area i, area j, kernel, state, align, session, prepost), each cell is a n_area_i x n_area_j matrix
J_data_err = cell(3, 3, n_conn_kernel, n_states, 3, n_session, 3);
firing_rates = cell(3, n_states, 3, n_session); % (area, state, prepost, session)

all_J = [];
J_state = [];
J_session = [];
J_kernel = [];
J_i = [];
J_j = [];

prepost_all = {'Pre', 'Post', 'Pre'}; % pre cortex, post cortex, pre full
for prepost_idx = 1:2
    prepost = prepost_all{prepost_idx};
    for session_idx = 1:n_session
        fprintf('Loading session %d\n', session_idx);
        for state_idx = 1:n_states
            state = states{state_idx};
            switch prepost_idx
                case 1
                    session_stage = [session_type, 'Pre', state, '_cortex'];
                case 2
                    session_stage = [session_type, 'Post', state, '_cortex'];
                case 3
                    session_stage = [session_type, 'Pre', state, '_full'];
            end

            for align_idx = 2:2
                align = aligns{align_idx};
                session_stage_full = [session_stage, '_', align];
                % file_path = [root_path, 'GLM_model/', session_stage_full,...
                %     '/GLM_', session_stage_full, '_', num2str(session_idx), '_',...
                %     kernel, '_0_', reg, '_', epoch, '.mat'];
        
                % fprintf('Loading %s\n', session_stage_full);
        
                % load(file_path, "model_par", "n_PS_kernel", "kernel_len", "N", "model_err");
                load([root_path, 'GLM_data/', session_stage_full,'/borders_', session_stage_full, '_', ...
                        num2str(session_idx),'.mat'], "borders");
                borders = [1, borders+0.5]; % area i is from borders(i) to borders(i+1)
                n_area = length(borders) - 1;

                file_path = [root_path, 'GLM_data/', session_stage_full,'/GLMdata_', session_stage_full, '_', ...
                        num2str(session_idx),'_DeltaPure_0.mat'];
                load(file_path, "raster");
                firing_rate = mean(raster, 2)*1000;
                for area_idx = 1:2
                    firing_rates{area_idx, state_idx, prepost_idx, session_idx} = firing_rate(borders(area_idx):borders(area_idx+1)-1);
                end
                % if strcmp(prepost, 'Post')
                %     assert(n_area == 2, 'Only 2 areas in Post sessions');
                % else
                %     assert(n_area == 3, 'Only 3 areas in Pre sessions');
                % end
            end
        end
    end
end

% 3x3 figure (area x state)
f = figure("Visible", "off", 'PaperPosition', [0 0 12 8]);
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% firing rate pre vs post
for state_idx = 1:n_states
    for area_idx = 1:2
        nexttile((area_idx-1)*3 + state_idx);
        % merge all sessions
        firing_rate_pre = [];
        firing_rate_post = [];
        for session_idx = 1:n_session
            firing_rate_pre = [firing_rate_pre; firing_rates{area_idx, state_idx, 1, session_idx}];
            firing_rate_post = [firing_rate_post; firing_rates{area_idx, state_idx, 2, session_idx}];
        end

        scatter(firing_rate_pre, firing_rate_post, 5, 'x');
        hold on;
        plot([0, 50], [0, 50], 'k--', 'LineWidth', 1);  % Diagonal line
        axis equal;
        xlim([0, 50]);
        ylim([0, 50]);
        set(gca, 'XTick', [0 25 50], 'YTick', [0 25 50]);
        xlabel('Pre');
        ylabel('Post');
        title([area_names{area_idx}, ' ', state_names{state_idx}]);

        % Add gray lines at zero
        line([0, 0], [0, 50], 'Color', 'black', 'LineWidth', 1);
        line([0, 50], [0, 0], 'Color', 'black', 'LineWidth', 1);
        hold off;
    end



end

% save figure
folder = [root_path, 'figures/firing_rate'];
check_path(folder);
saveas(f, [folder, '/firing_rate_', session_type, '.png']);
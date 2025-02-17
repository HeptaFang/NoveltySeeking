session_type = 'Muscimol';
root_path = '../';
kernel = 'DeltaPure';
reg = 'L2=2';
epoch = '2500';
area_names = {'ACC', 'Thalamus', 'VLPFC'};
states = {'Task', 'RestOpen', 'RestClose'};
aligns = {'AlignLast'};
n_states = length(states);
significant_threshold = 1;
use_abs_mean = false;
use_abs_sigmean = false;
% sign_filter = 'raw'; % 'abs', 'pos', 'neg', or 'raw'
plot_legend = true;

if strcmp(session_type, 'Muscimol')
    n_session = 10;
else
    n_session = 5;
end
n_conn_kernel = 3;

% warning('off');


prepost_all = {'Pre', 'Post'};
for prepost_idx = 1:2
    prepost = prepost_all{prepost_idx};
    for session_idx = 1:n_session
        fprintf('Loading session %d\n', session_idx);
        for state_idx = 1:n_states
            state = states{state_idx};
            if strcmp(prepost, 'Pre')
                session_stage = [session_type, 'Pre', state, '_full'];
                % session_stage_full = [session_type, 'Pre', state, '_cortex'];
            else
                session_stage = [session_type, 'Post', state, '_cortex'];
            end

            for align_idx = 1:1
                align = aligns{align_idx};
                session_stage_full = [session_stage, '_', align];
                file_path = [root_path, 'GLM_model/', session_stage_full,...
                    '/GLM_', session_stage_full, '_', num2str(session_idx), '_',...
                    kernel, '_0_', reg, '_', epoch, '.mat'];
        
                fprintf('Loading %s\n', session_stage_full);
        
                load(file_path, "model_par", "n_PS_kernel", "kernel_len", "N", "model_err");
                load([root_path, 'GLM_data/', session_stage_full,'/borders_', session_stage_full, '_', ...
                        num2str(session_idx),'.mat'], "borders");

                file_path = [root_path, 'GLM_data/', session_stage_full,...
                    '/GLMdata_', session_stage_full, '_', num2str(session_idx), '_',...
                    kernel, '_0.mat'];
                load(file_path, "raster")
                firing_rate = mean(raster, 2) * 1000;
                firing_rates{state_idx, align_idx, session_idx, prepost_idx} = firing_rate;
                
                borders = [1, borders+0.5]; % area i is from borders(i) to borders(i+1)
                n_area = length(borders) - 1;
                % if strcmp(prepost, 'Post')
                %     assert(n_area == 2, 'Only 2 areas in Post sessions');
                % else
                %     assert(n_area == 3, 'Only 3 areas in Pre sessions');
                % end
                f = figure('Position', [100, 100, 1200, 800], "Visible", "off");
                t = tiledlayout(2, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
        
                for kernel_idx = 1:n_conn_kernel
                    ax = nexttile(t, kernel_idx);
                    J_mat = model_par(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));
                    J_mean_in = mean(J_mat, 1);
                    J_mean_out = mean(J_mat, 2);
                    hold on;
                    scatter(firing_rate, J_mean_in, 10, 'x');
                    scatter(firing_rate, J_mean_out, 10, '+');
                    hold off;
                    xlabel('Firing rate (Hz)');
                    ylabel('average J');
                    title(['Kernel', num2str(kernel_idx)]);
                    legend('J in', 'J out');
                end
                for kernel_idx = 1:n_conn_kernel
                    ax = nexttile(t, kernel_idx+3);
                    J_mat = model_par(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));
                    J_mat = abs(J_mat);
                    J_mean_in = mean(J_mat, 1);
                    J_mean_out = mean(J_mat, 2);
                    hold on;
                    scatter(firing_rate, J_mean_in, 10, 'x');
                    scatter(firing_rate, J_mean_out, 10, '+');
                    hold off;
                    xlabel('Firing rate (Hz)');
                    ylabel('average abs J');
                    title(['Kernel', num2str(kernel_idx)]);
                    legend('J in', 'J out');
                end



                title(t, [session_stage_full, ' ', num2str(session_idx)]);
                save_folder = [root_path, 'figures/fr_conn/', session_stage_full]; 
                check_path(save_folder);
                saveas(f, [save_folder, '/fr_conn_', session_stage_full, '_', num2str(session_idx), '_',...
                    kernel,'.png']);
            end
        end
    end
end
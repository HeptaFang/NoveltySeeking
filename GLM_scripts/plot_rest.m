function plot_rest(session, kernel_name, epoch, reg, shuffle_size)
%% open vs close only, for each connection
datasets = {'eyeOpen', 'eyeClose'};
open_data = struct();
close_data = struct();
for data_idx=1:2
    dataset_name = datasets{data_idx};
    % load model parameters
    model_path_ori = ['../GLM_model/', dataset_name, '/GLM_', dataset_name, '_', ...
            int2str(session), '_', kernel_name, '_0_', ...
            reg.name, '_', int2str(epoch)];
    load(model_path_ori, "model_par", "PS_kernels", "conn_kernels", "n_PS_kernel", "n_conn_kernel", "kernel_len", "N");
    
    if session<10
        session_border = session;
    else
        session_border = floor(session/100);
    end
    
    load(['../GLM_data/', dataset_name,'/borders_', dataset_name, '_', ...
            int2str(session_border),'.mat'], "borders");
    % A_T_border = borders(1);
    % T_P_border = borders(2);
    borders = [0, borders-0.5, N];
    
    par_ori = model_par;
    par_sfl = zeros([size(par_ori), shuffle_size]);
    
    for i=1:shuffle_size
        model_path_sfl = ['../GLM_model/', dataset_name, '/GLM_', dataset_name, '_', ...
            int2str(session), '_', kernel_name, '_', int2str(i), '_', ...
            reg.name, '_', int2str(epoch)];
        par_sfl(:, :, i) = load(model_path_sfl).model_par;
    end

    % extract J of 1st kernel
    par_ori = par_ori(:, 2:end);
    par_sfl = par_sfl(:, 2:end, :);

    % remove cell #101 in session 4 eye_close. (#45 in aligned dataset)
    if session==4 && data_idx==2
        par_ori = [par_ori(:, 1:44),par_ori(:, 46:end)];
        par_ori = [par_ori(1:44, :);par_ori(46:end, :)];
        par_sfl = [par_sfl(:, 1:44, :),par_sfl(:, 46:end, :)];
        par_sfl = [par_sfl(1:44, :, :);par_sfl(46:end, :, :)];
        borders(end) = borders(end)-1;
    end
    
    % split blocks
    J_block_ori = cell(3,3);
    J_block_sfl = cell(3,3);
    
    for i = 1:3
        for j = 1:3
            i_range = (borders(i)+1):(borders(i+1));
            j_range = (borders(j)+1):(borders(j+1));
            J_block_ori{i, j} = par_ori(i_range, j_range);
            J_block_sfl{i, j} = par_sfl(i_range, j_range, :);
        end
    end
    if data_idx==1
        open_data.ori=J_block_ori;
        open_data.sfl=J_block_sfl;
    else
        close_data.ori=J_block_ori;
        close_data.sfl=J_block_sfl;
    end
end


% plot
fig = figure("Visible", "off");
set(fig, 'PaperPosition', [0, 0, 30, 30]);
tiles = tiledlayout(3, 6);
cmap="jet";

% y: kernel, ori, significant, shuffle1, ..., shuffleN
% x: h, conn_kernels, PS_kernels
areas = {'A', 'T', 'V'};

for plot_x=1:6
    for plot_y = 1:3
        t = nexttile((plot_y-1)*6+plot_x);

        if mod(plot_x, 2)==1
            % original plot
            open = open_data.ori;
            close = close_data.ori;
            
        else
            % shuffle (ave over shuffles?)
            open = open_data.sfl;
            close = close_data.sfl;
        end
        
        open = open{plot_y, ceil(plot_x/2)};
        close = close{plot_y, ceil(plot_x/2)};
        open = reshape(open, numel(open), 1);
        close = reshape(close, numel(close), 1);

        if length(open)~=length(close)
            fprintf("length mismatch\n");
            continue
        end

        % plot
        plot([open, close].', '.-', 'Color', [1, 0.6, 0.6]);
        xticks([1, 2]);
        xticklabels({'Open', 'Close'});
        xlim([0.5,2.5]);
        ylim([-2, 2]);

        if length(open)>0
            ave = [mean(open), mean(close)];
            err = [std(open), std(close)];
            hold on;
            errorbar(ave, err, "Color", [0, 0, 0],"LineWidth",1.5);

        if mod(plot_x, 2)==1
            title(sprintf("%s-%s ori", areas{ceil(plot_x/2)}, areas{plot_y}));
        else
            title(sprintf("%s-%s sfl", areas{ceil(plot_x/2)}, areas{plot_y}));
        end

    end
end

fig_path = ['../figures/GLM/', 'restStates'];
check_path(fig_path);
fig_file = [fig_path, '/restStates_' ...
        int2str(session), '_', kernel_name, '_', ...
        reg.name, '_', int2str(epoch), '.png'];
% exportgraphics(fig,fig_file);
print(fig, fig_file,'-dpng', '-r300');
end
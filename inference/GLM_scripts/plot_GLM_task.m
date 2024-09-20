tasks = {'cue_shuffle'};

% model_suffix = '_shuffleAll';
model_suffix = '';
data = load("combined_aligned.mat").aligned_data;
align_idxs = {data.spk_open_idx; data.spk_close_idx};
for task_idx = 1:1
    align_idx_cell = align_idxs{task_idx};
    align_idx = cellfun(@mean, align_idx_cell);
    for session_idx = 1:5
        for epoch = 100:100:1000
            taskname=[tasks{task_idx}];
            session_name=['session', int2str(session_idx)];
            filename = ['GLM_results/',taskname,'/',session_name,'_',int2str(epoch),model_suffix];
            load(filename); % par_vector, minuslogL, N
            
            h = par_vector(1:N);
            J = zeros(N, N);
            
            for i=1:N
                J(i, 1:i-1) = par_vector(N + (i-1)*(N-1)+(1:(i-1)), 1);
                J(i, i+1:N) = par_vector(N + (i-1)*(N-1)+(i:(N-1)), 1);
            end
            
            fig = figure("Visible","off");
            subplot(1, 2, 1);
            imagesc(J);
            axis square
            colormap(jet);
            colorbar;
            xlabel('From cell idx');
            ylabel('To cell idx');
            clim([-3, 3]);
            title(sprintf("J, %s, session%d", tasks{task_idx}, session_idx));
    
            A_T_border = find(align_idx(:, session_idx)>45, 1);
            T_P_border = find(align_idx(:, session_idx)>98, 1);
            if isempty(A_T_border)
                A_T_border = N + 1;
            end
            if isempty(T_P_border)
                T_P_border = N + 1;
            end
            A_T_border = A_T_border-0.5;
            T_P_border = T_P_border-0.5;
    
            hold on
            if A_T_border == T_P_border
                line([0,N+0.5], [A_T_border,A_T_border], 'Color', 'k');
                line([A_T_border,A_T_border], [0,N+0.5], 'Color', 'k');
            else
                line([0,N+0.5], [A_T_border,A_T_border], 'Color', 'r');
                line([0,N+0.5], [T_P_border,T_P_border], 'Color', 'b');
                line([A_T_border,A_T_border], [0,N+0.5], 'Color', 'r');
                line([T_P_border,T_P_border], [0,N+0.5], 'Color', 'b');
            end
            hold off
            save(sprintf("GLMresult_%s_%d_%d.mat", tasks{task_idx}, session_idx, epoch), "J", "h", "A_T_border", "T_P_border");
            
            subplot(1, 2, 2);
            plot(1:N, h);
            xlabel('Cell idx');
            ylabel('h');
            title(sprintf("h"));
            
            folderName = ['GLM_results/',taskname];
            if ~exist(folderName, 'dir')
                mkdir(folderName);
            end
            saveas(fig, sprintf("figures/GLM_directed_%s_%d%s_%d.png", tasks{task_idx}, session_idx, model_suffix, epoch));
            close(fig);
        end
    end
end


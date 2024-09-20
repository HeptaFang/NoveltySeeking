
ACC_data = load("ACC.mat").ACC;
thalamus_data = load("thalamus.mat").thalamus;
VLPFC_data = load("VLPFC.mat").VLPFC;
% combined_data = load("combined.mat");
% ACC_data = load("ACC_aligned.mat").aligned_data;
% thalamus_data = load("thalamus_aligned.mat").aligned_data;
% VLPFC_data = load("VLPFC_aligned.mat").aligned_data;
combined_data = load("combined_aligned.mat").aligned_data;

combined_uncID = load("combined_uncID.mat").combined_uncID;

datas = {ACC_data; thalamus_data; VLPFC_data; combined_data};
data_names = {"ACC"; "thalamus"; "VLPFC"; "combined"};

% List of inference output filenames
filenames = {'ACC_eye_open', 'ACC_eye_close','thalamus_eye_open', 'thalamus_eye_close','VLPFC_eye_open', 'VLPFC_eye_close', 'combined_eye_open', 'combined_eye_close',}; 
% filenames = {'combined_eye_open',}; 

% Loop over each file
for data_idx = 4:length(datas)
    data = datas{data_idx};
    disp(data_names{data_idx});

    subdatasets = {data.spk_open; data.spk_close; data.spk_iti; data.spk_cue};
    align_idxs = {data.spk_open_idx; data.spk_close_idx;};
    sub_names = {"eye open"; "eye close"; "ITI"; "cue"};

    ave_J = zeros(5, 3, 2, 2);
    std_J = zeros(5, 3, 2, 2);

    % no task
    for session_idx = 1:5
        % f = figure('visible','off');
        fig = figure('visible', 'off');
        skipped = false;
        open_data = [];
        close_data = [];
        for sub_idx = 1:2
            subdataset = subdatasets{sub_idx};
            align_idx_cell = align_idxs{sub_idx};
            align_idx = cellfun(@mean, align_idx_cell);
            original_data = subdataset(:, session_idx);
            valid = ones(1, length(original_data));
            max_n = 0;
            for i = 1:length(original_data)
                if isempty(original_data{i})
                    valid(i) = 0;
                else
                    max_n = i;
                end
            end
            
            fprintf("session %d, valid:%d\n", session_idx, max_n)
            if max_n == 0
                skipped = true;
                continue
            end

            fileIdx = (data_idx - 1) * 2 + sub_idx;
            
            plots = ["_bin5.p", "_bin5.p", "_L2_bin5.j"];

            file_path=['./', 'inference_input/', filenames{fileIdx}, '/saved_data/parameters/overlapping/'];
            filename_moment = [file_path, 'session', int2str(session_idx), '_bin5.p'];
            filename_connection = [file_path, 'session', int2str(session_idx), '_L2_bin5.j'];
            filename_error = [file_path, 'session', int2str(session_idx), '_L2_bin5_err.j'];
        
            % Read the file
            skip_flag=false;
            
            fid_moment = fopen(filename_moment, 'r');
            fid_connection = fopen(filename_connection, 'r');
            fid_error = fopen(filename_error, 'r');
            moment_exist = (fid_error == -1);
            connection_exist = (fid_error == -1);
            error_exist = (fid_error ~= -1);
            fprintf("Error: %d\n", error_exist);

            if fid_connection == -1
                fprintf("data %s %d skip\n", filenames{fileIdx}, session_idx)
                skip_flag=true;
            end
            fprintf("data %s session%d start\n", filenames{fileIdx}, session_idx)
            moment_data = fscanf(fid_moment, '%f');
            fclose(fid_moment);
            if ~skip_flag
                connection_data = fscanf(fid_connection, '%f');
                fclose(fid_connection);
            end
            if error_exist
                error_data = importdata(filename_error);
                error_data = error_data(:, 2);
            end
        
            % Determine N
            numLines = length(moment_data);
            if ~skip_flag
                numLines_c = length(connection_data);
                
                if numLines_c ~= numLines
                    fprintf("line number not match: moment:%d, connection:%d.", numLines. numLines_c);
                end
            end
            N = (1 + sqrt(1 + 8 * numLines)) / 2;
            if floor(N) ~= N
                error('Invalid format in file %s', filename);
            end
            disp(N)
            disp(max_n)
            
            % only for unaligned
            max_n = N-1;
        
            % Construct the connection matrix
            corr = ones(max_n, max_n) * NaN;
            J = ones(max_n, max_n) * NaN;
            err = ones(max_n, max_n) * NaN;
            corr_flat = zeros(size(moment_data));
            index = 1;
            for i = 1:max_n
                for j = i+1:max_n
                    if valid(i) && valid(j)
                        corr(i, j) = moment_data(index+max_n)/(moment_data(i)*moment_data(j));
                        corr(j, i) = moment_data(index+max_n)/(moment_data(i)*moment_data(j));
                        corr_flat(index+max_n) = moment_data(index+max_n)/(moment_data(i)*moment_data(j));
                        if ~skip_flag
                            J(i, j) = connection_data(index+max_n);
                            J(j, i) = connection_data(index+max_n);
                        end
                        if error_exist
                            err(i, j) = error_data(index+max_n);
                            err(j, i) = error_data(index+max_n);
                        end
                    else
                        corr(i, j) = NaN;
                        corr(j, i) = NaN;
                        if ~skip_flag
                            J(i, j) = NaN;
                            J(j, i) = NaN;
                        end
                    end
                    index = index + 1;
                end
            end

            if sub_idx == 1
                open_data = J;
                if error_exist
                    open_err = err;
                end
            end
            if sub_idx == 2 
                close_data = J;
                if error_exist
                    close_err = err;
                end
            end
             

            % Plot
            subplot(3, 3, (sub_idx-1)*3+1);
            h = imagesc(log(corr));
            axis square
            set(h, 'AlphaData', ~isnan(corr))
            colormap(jet);
            colorbar;
            xlabel('Cell idx');
            ylabel('Cell idx');
            title(sprintf("ln(Corr idx), %s", sub_names{sub_idx}));
            clim([-2, 2])
            if data_idx == 4
                A_T_border = find(align_idx(:, session_idx)>45, 1);
                T_P_border = find(align_idx(:, session_idx)>98, 1);
                if isempty(A_T_border)
                    A_T_border = max_n + 1;
                end
                if isempty(T_P_border)
                    T_P_border = max_n + 1;
                end
                A_T_border = A_T_border-0.5;
                T_P_border = T_P_border-0.5;

                hold on
                if A_T_border == T_P_border
                    line([0,max_n+0.5], [A_T_border,A_T_border], 'Color', 'k');
                    line([A_T_border,A_T_border], [0,max_n+0.5], 'Color', 'k');
                else
                    line([0,max_n+0.5], [A_T_border,A_T_border], 'Color', 'r');
                    line([0,max_n+0.5], [T_P_border,T_P_border], 'Color', 'b');
                    line([A_T_border,A_T_border], [0,max_n+0.5], 'Color', 'r');
                    line([T_P_border,T_P_border], [0,max_n+0.5], 'Color', 'b');
                end
                hold off
            end
            
            if ~skip_flag
                subplot(3, 3, (sub_idx-1)*3+2);
                h = imagesc(J);
                axis square
                set(h, 'AlphaData', ~isnan(J))
                colormap(jet);
                colorbar;
                xlabel('Cell idx');
                ylabel('Cell idx');
                title(sprintf("J_{ij}, %s", sub_names{sub_idx}));
                clim([-2, 2])
                if data_idx == 4
                    hold on

                    if A_T_border == T_P_border
                        line([0,max_n+0.5], [A_T_border,A_T_border], 'Color', 'k');
                        line([A_T_border,A_T_border], [0,max_n+0.5], 'Color', 'k');
                    else
                        line([0,max_n+0.5], [A_T_border,A_T_border], 'Color', 'r');
                        line([0,max_n+0.5], [T_P_border,T_P_border], 'Color', 'b');
                        line([A_T_border,A_T_border], [0,max_n+0.5], 'Color', 'r');
                        line([T_P_border,T_P_border], [0,max_n+0.5], 'Color', 'b');
                    end
                    hold off
                end

                subplot(3, 3, (sub_idx-1)*3+3);
                h = imagesc(J);
                axis square
                set(h, 'AlphaData', ~isnan(J))
                colormap(jet);
                colorbar;
                xlabel('Cell idx');
                ylabel('Cell idx');
                title(sprintf("significant J_{ij}, %s", sub_names{sub_idx}));
                clim([-2, 2])

                % mask
                hold on
                if error_exist
                    mask = abs(J)>2*abs(err);
                    [mask_x, mask_y] = find(mask);
                    scatter(mask_x, mask_y, 8, 'black');
                end
                hold off

                if data_idx == 4
                    hold on

                    if A_T_border == T_P_border
                        line([0,max_n+0.5], [A_T_border,A_T_border], 'Color', 'k');
                        line([A_T_border,A_T_border], [0,max_n+0.5], 'Color', 'k');
                    else
                        line([0,max_n+0.5], [A_T_border,A_T_border], 'Color', 'r');
                        line([0,max_n+0.5], [T_P_border,T_P_border], 'Color', 'b');
                        line([A_T_border,A_T_border], [0,max_n+0.5], 'Color', 'r');
                        line([T_P_border,T_P_border], [0,max_n+0.5], 'Color', 'b');
                    end
                    hold off
                end

                subplot(3, 3, sub_idx+6);
                scatter(log(corr_flat(max_n+1:end)), connection_data(max_n+1:end));
                all_data = [log(corr_flat(max_n+1:end)), connection_data(max_n+1:end)];
                lim = [min(all_data(~isinf(all_data)), [], 'all'), ...
                    max(all_data(~isinf(all_data)), [], 'all')];
                disp(lim)
                hold on
                xy_line = plot(lim, lim);
                legend(xy_line, 'x=y', "Location","bestoutside");
                hold off
                xlabel('ln(Corr idx)');
                ylabel('J');
                title(sprintf("%s", sub_names{sub_idx}));
            end

        end
        if ~skipped
            % Save the figure
            % todo: plot diff.

            subplot(3, 3, 9);
            if data_idx == 4
                % remove cell #101 in session 4 eye_close. (#45 in aligned dataset)
                if session_idx == 4
                    close_data = [close_data(:, 1:44),close_data(:, 46:end)];
                    close_data = [close_data(1:44, :);close_data(46:end, :)];
                    close_err = [close_err(:, 1:44),close_err(:, 46:end)];
                    close_err = [close_err(1:44, :);close_err(46:end, :)];
                    max_n = max_n - 1;
                end

                conn_group = zeros(max_n);
                mask = (abs(open_data)>2*abs(open_err)) + (abs(close_data)>2*abs(close_err)) > 0;
                align_idx_valid = align_idx(:, session_idx);
                align_idx_valid = align_idx_valid(~isnan(align_idx_valid));
                uncID = combined_uncID{1,session_idx}(align_idx_valid);
                unc_group = zeros(max_n);
                align_idx;

                for i = 1:max_n
                    for j = 1:max_n
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
                        end
                        if i < A_T_border && j >= T_P_border
                            conn_group(i, j) = 5; % ACC-VLPFC, M
                        end
                        if i >= A_T_border && i < T_P_border  && j >= T_P_border
                            conn_group(i, j) = 6; % Thalamus-VLPFC, C
                        end
                        
                        unc_group(i, j) = uncID(i) + uncID(j);

                    end
                end
                % plot
                maxlim = max([open_data, close_data], [], 'all');
                minlim = min([open_data, close_data], [], 'all');
                lim = [minlim, maxlim];
                xy_line = plot(lim, lim);
                hold on

                % idx = find(((conn_group == 1) + ~mask)==2);
                % scatter(open_data(idx), close_data(idx), "red", ".");
                % idx = find(((conn_group == 2) + ~mask)==2);
                % scatter(open_data(idx), close_data(idx), "green", ".");
                % idx = find(((conn_group == 3) + ~mask)==2);
                % scatter(open_data(idx), close_data(idx), "blue", ".");
                % idx = find(((conn_group == 4) + ~mask)==2);
                % scatter(open_data(idx), close_data(idx), "yellow", ".");
                % idx = find(((conn_group == 5) + ~mask)==2);
                % scatter(open_data(idx), close_data(idx), "magenta", ".");
                % idx = find(((conn_group == 6) + ~mask)==2);
                % scatter(open_data(idx), close_data(idx), "cyan", ".");
                
                mkr = '.';
                % idx = find(((conn_group == 1) + mask)==2);
                % scatter(open_data(idx), close_data(idx), "red", mkr);
                % idx = find(((conn_group == 2) + mask)==2);
                % scatter(open_data(idx), close_data(idx), "green", mkr);
                % idx = find(((conn_group == 3) + mask)==2);
                % scatter(open_data(idx), close_data(idx), "blue", mkr);
                % idx = find(((conn_group == 4) + mask)==2);
                % scatter(open_data(idx), close_data(idx), "yellow", mkr);
                % idx = find(((conn_group == 5) + mask)==2);
                % scatter(open_data(idx), close_data(idx), "magenta", mkr);
                % idx = find(((conn_group == 6) + mask)==2);
                % scatter(open_data(idx), close_data(idx), "cyan", mkr);

                idx = find(unc_group==0);
                scatter(open_data(idx), close_data(idx), "red", mkr);
                idx = find(unc_group==1);
                scatter(open_data(idx), close_data(idx), "magenta", mkr);
                idx = find(unc_group==2);
                scatter(open_data(idx), close_data(idx), "blue", mkr);
                hold off

                % legend(xy_line, 'x=y', "Location","bestoutside");
                % legend(["x=y", "ACC-ACC", "Tha-Tha", "PFC-PFC", "ACC-Tha", "ACC-PFC", "Tha-PFC"], "Location","bestoutside");
                % legend(["x=y", "A-A", "T-T", "P-P", "A-T", "A-P", "T-P"], "Location","bestoutside");
                legend(["x=y", "N-N", "U-N", "U-U"], "Location","bestoutside");
                % title('significant J');
                title('J');
                xlabel('J_{ij} eye open');
                ylabel('J_{ij} eye close');
                set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
            end
            
            sgtitle(sprintf("%s, session%d", data_names{data_idx}, session_idx))
            % print(fig, sprintf('figures/%s_session%d.png', data_names{data_idx}, session_idx), 'Resolution', 600);
            fig.Position = [0 0 1200 1000];
            saveas(fig, sprintf('figures/%s_session%d.png', data_names{data_idx}, session_idx));

            % single Hi-res fig
            fig = figure('Visible', 'off');
            maxlim = max([open_data, close_data], [], 'all');
            minlim = min([open_data, close_data], [], 'all');
            lim = [minlim, maxlim];
            xy_line = plot(lim, lim);
            hold on

            mkr = '.';
            % idx = find(((conn_group == 1) + ~mask)==2);
            % scatter(open_data(idx), close_data(idx), "red", mkr);
            % idx = find(((conn_group == 2) + ~mask)==2);
            % scatter(open_data(idx), close_data(idx), "green", mkr);
            % idx = find(((conn_group == 3) + ~mask)==2);
            % scatter(open_data(idx), close_data(idx), "blue", mkr);
            % idx = find(((conn_group == 4) + ~mask)==2);
            % scatter(open_data(idx), close_data(idx), "yellow", mkr);
            % idx = find(((conn_group == 5) + ~mask)==2);
            % scatter(open_data(idx), close_data(idx), "magenta", mkr);
            % idx = find(((conn_group == 6) + ~mask)==2);
            % scatter(open_data(idx), close_data(idx), "cyan", mkr);
            idx = find((unc_group==0) & (~mask));
            scatter(open_data(idx), close_data(idx), "red", mkr);
            idx = find((unc_group==1) & (~mask));
            scatter(open_data(idx), close_data(idx), "magenta", mkr);
            idx = find((unc_group==2) & (~mask));
            scatter(open_data(idx), close_data(idx), "blue", mkr);
            
            mkr = 'o';
            % idx = find(((conn_group == 1) + mask)==2);
            % scatter(open_data(idx), close_data(idx), "red", mkr);
            % idx = find(((conn_group == 2) + mask)==2);
            % scatter(open_data(idx), close_data(idx), "green", mkr);
            % idx = find(((conn_group == 3) + mask)==2);
            % scatter(open_data(idx), close_data(idx), "blue", mkr);
            % idx = find(((conn_group == 4) + mask)==2);
            % scatter(open_data(idx), close_data(idx), "yellow", mkr);
            % idx = find(((conn_group == 5) + mask)==2);
            % scatter(open_data(idx), close_data(idx), "magenta", mkr);
            % idx = find(((conn_group == 6) + mask)==2);
            % scatter(open_data(idx), close_data(idx), "cyan", mkr);
            idx = find((unc_group==0) & (mask));
            scatter(open_data(idx), close_data(idx), "red", mkr);
            idx = find((unc_group==1) & (mask));
            scatter(open_data(idx), close_data(idx), "magenta", mkr);
            idx = find((unc_group==2) & (mask));
            scatter(open_data(idx), close_data(idx), "blue", mkr);

            %%%%%%% specfic area
            mask = mask & (conn_group==3);
            %%%%%%%

            idx = find((unc_group==0) & (mask) & open_data > 0);
            ave_J(session_idx, 1, 1, 1) = mean(open_data(idx));
            std_J(session_idx, 1, 1, 1) = std(open_data(idx));
            idx = find((unc_group==0) & (mask) & close_data > 0);
            ave_J(session_idx, 1, 2, 1) = mean(close_data(idx));
            std_J(session_idx, 1, 2, 1) = std(close_data(idx));
            idx = find((unc_group==0) & (mask) & open_data < 0);
            ave_J(session_idx, 1, 1, 2) = mean(open_data(idx));
            std_J(session_idx, 1, 1, 2) = std(open_data(idx));
            idx = find((unc_group==0) & (mask) & close_data < 0);
            ave_J(session_idx, 1, 2, 2) = mean(close_data(idx));
            std_J(session_idx, 1, 2, 2) = std(close_data(idx));

            idx = find((unc_group==1) & (mask) & open_data > 0);
            ave_J(session_idx, 2, 1, 1) = mean(open_data(idx));
            std_J(session_idx, 2, 1, 1) = std(open_data(idx));
            idx = find((unc_group==1) & (mask) & close_data > 0);
            ave_J(session_idx, 2, 2, 1) = mean(close_data(idx));
            std_J(session_idx, 2, 2, 1) = std(close_data(idx));
            idx = find((unc_group==1) & (mask) & open_data < 0);
            ave_J(session_idx, 2, 1, 2) = mean(open_data(idx));
            std_J(session_idx, 2, 1, 2) = std(open_data(idx));
            idx = find((unc_group==1) & (mask) & close_data < 0);
            ave_J(session_idx, 2, 2, 2) = mean(close_data(idx));
            std_J(session_idx, 2, 2, 2) = std(close_data(idx));
            
            idx = find((unc_group==2) & (mask) & open_data > 0);
            ave_J(session_idx, 3, 1, 1) = mean(open_data(idx));
            std_J(session_idx, 3, 1, 1) = std(open_data(idx));
            idx = find((unc_group==2) & (mask) & close_data > 0);
            ave_J(session_idx, 3, 2, 1) = mean(close_data(idx));
            std_J(session_idx, 3, 2, 1) = std(close_data(idx));
            idx = find((unc_group==2) & (mask) & open_data < 0);
            ave_J(session_idx, 3, 1, 2) = mean(open_data(idx));
            std_J(session_idx, 3, 1, 2) = std(open_data(idx));
            idx = find((unc_group==2) & (mask) & close_data < 0);
            ave_J(session_idx, 3, 2, 2) = mean(close_data(idx));
            std_J(session_idx, 3, 2, 2) = std(close_data(idx));

            ax = gca;
            set(ax, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
            
            set(gcf, 'InvertHardCopy', 'off'); 
            set(ax, 'color', [0.8 0.8 0.8]);
            set(gcf, 'color', [0.9 0.9 0.9]);

            hold off
            % legend(xy_line, 'x=y', "Location","bestoutside");
            % legend(["x=y", "ACC-ACC", "Tha-Tha", "PFC-PFC", "ACC-Tha", "ACC-PFC", "Tha-PFC"], "Location","bestoutside");
            % legend(["x=y", "A-A", "T-T", "P-P", "A-T", "A-P", "T-P"], "Location","bestoutside");
            legend(["x=y", "NonUnc-NonUnc", "Unc-NonUnc", "Unc-Unc"], "Location","bestoutside");
            title('J');
            xlabel('J_{ij} eye open');
            ylabel('J_{ij} eye close');
            fig.Position = [0 0 1200 1000];
            
            saveas(fig, sprintf('figures/%s_session%d_diff.png', data_names{data_idx}, session_idx));

        end
        close all;
    end

    fig = figure('Visible', 'off');
    
    jit_size=0.25;
    
    errorbar((1:5)-0.5*jit_size, ave_J(:, 1, 1, 1), std_J(:, 1, 1, 1), "Marker","+",LineStyle="-", Color='r');
    hold on;
    errorbar((1:5)-0.3*jit_size, ave_J(:, 2, 1, 1), std_J(:, 2, 1, 1), "Marker","x",LineStyle="-", Color='m');
    errorbar((1:5)-0.1*jit_size, ave_J(:, 3, 1, 1), std_J(:, 3, 1, 1), "Marker","o",LineStyle="-", Color='b');
    errorbar((1:5)+0.1*jit_size, ave_J(:, 1, 2, 1), std_J(:, 1, 2, 1), "Marker","+",LineStyle="--", Color='r');
    errorbar((1:5)+0.3*jit_size, ave_J(:, 2, 2, 1), std_J(:, 2, 2, 1), "Marker","x",LineStyle="--", Color='m');
    errorbar((1:5)+0.5*jit_size, ave_J(:, 3, 2, 1), std_J(:, 3, 2, 1), "Marker","o",LineStyle="--", Color='b');

    errorbar((1:5)-0.5*jit_size, ave_J(:, 1, 1, 2), std_J(:, 1, 2, 2), "Marker","+",LineStyle="-", Color='r');
    errorbar((1:5)-0.3*jit_size, ave_J(:, 2, 1, 2), std_J(:, 2, 1, 2), "Marker","x",LineStyle="-", Color='m');
    errorbar((1:5)-0.1*jit_size, ave_J(:, 3, 1, 2), std_J(:, 3, 1, 2), "Marker","o",LineStyle="-", Color='b');
    errorbar((1:5)+0.1*jit_size, ave_J(:, 1, 2, 2), std_J(:, 1, 2, 2), "Marker","+",LineStyle="--", Color='r');
    errorbar((1:5)+0.3*jit_size, ave_J(:, 2, 2, 2), std_J(:, 2, 2, 2), "Marker","x",LineStyle="--", Color='m');
    errorbar((1:5)+0.5*jit_size, ave_J(:, 3, 2, 2), std_J(:, 3, 2, 2), "Marker","o",LineStyle="--", Color='b');
    hold off;

    fig.Position = [0 0 1200 1000];
    title('J VLPFC');
    xlabel('session');
    ylabel('J');
    legend("NonUnc-NonUnc open", "Unc-NonUnc open", "Unc-Unc open", ...
        "NonUnc-NonUnc close", "Unc-NonUnc close", "Unc-Unc close");
    saveas(fig, sprintf('figures/%s_UncEnc_err_VLPFC.png', data_names{data_idx}));

    fig = figure('Visible', 'off');
    plot((1:5)-0.5*jit_size, ave_J(:, 1, 1), "Marker","+",LineStyle="-", Color='r');
    hold on;
    plot((1:5)-0.3*jit_size, ave_J(:, 2, 1), "Marker","x",LineStyle="-", Color='m');
    plot((1:5)-0.1*jit_size, ave_J(:, 3, 1), "Marker","o",LineStyle="-", Color='b');
    plot((1:5)+0.1*jit_size, ave_J(:, 1, 2), "Marker","+",LineStyle="--", Color='r');
    plot((1:5)+0.3*jit_size, ave_J(:, 2, 2), "Marker","x",LineStyle="--", Color='m');
    plot((1:5)+0.5*jit_size, ave_J(:, 3, 2), "Marker","o",LineStyle="--", Color='b');
    hold off;

    legend("NonUnc-NonUnc open", "Unc-NonUnc open", "Unc-Unc open", ...
        "NonUnc-NonUnc close", "Unc-NonUnc close", "Unc-Unc close");

    fig.Position = [0 0 1200 1000];
    title('J');
    xlabel('session');
    ylabel('J');

    saveas(fig, sprintf('figures/%s_UncEnc.png', data_names{data_idx}));
end




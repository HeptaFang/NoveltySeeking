function kernel_weight_comparason(dataset_name, session, kernel_name, epoch, reg, shuffle_size, sorting, idx)
model_path_ori = ['../GLM_model/', dataset_name, '/GLM_', dataset_name, '_', ...
        int2str(session), '_', kernel_name, '_0_', ...
        reg.name, '_', int2str(epoch)];
load(model_path_ori, "model_par", "PS_kernels", "conn_kernels", "n_PS_kernel", "n_conn_kernel", "kernel_len", "N");

data_path = ['../GLM_data/', dataset_name, '/GLMdata_', dataset_name, '_', ...
        int2str(session),'_', kernel_name, '_0.mat'];
load(data_path, "raster");
firing_rate = mean(raster, 2);

if session<100
    session_border = session;
else
    session_border = floor(session/100);
end

load(['../GLM_data/', dataset_name,'/borders_', dataset_name, '_', ...
        int2str(session_border),'.mat'], "borders");

% for acc-thalamus-vlpfc dataset
A_T_border = borders(1);
T_P_border = borders(2);
    
par_ori = model_par;

par_sfl = zeros([size(par_ori), shuffle_size]);
for i=1:shuffle_size
    model_path_sfl = ['../GLM_model/', dataset_name, '/GLM_', dataset_name, '_', ...
        int2str(session), '_', kernel_name, '_', int2str(i), '_', ...
        reg.name, '_', int2str(epoch)];
    par_sfl(:, :, i) = load(model_path_sfl).model_par;
end

significant = zeros(size(par_ori));

for i=1:N
    for j=1:(1+n_PS_kernel+N*n_conn_kernel)
        % one-sample t-test
        h = ttest(reshape(par_sfl(i, j, :)-par_ori(i, j), [1, shuffle_size]), ...
            0, "Alpha", 0.001);
        significant(i, j) = h;
        % permutation test (doesn't work!)
        % p = permutationTest(par_ori(i, j), reshape(par_sfl(i, j, :), [1, shuffle_size]), 0, 'exact', 1);
        % significant(i, j) = p<0.05;

    end
end



par_sig = par_ori;
par_sig = (par_sig - mean(par_sfl, 3)).*significant;
save("par_sig.mat", "par_sig");
% par_sig(~significant) = nan;

% clim for J_ijk plottings
clim_ori = max(abs(par_ori(:, (n_PS_kernel+2):end)), [], "all");
clim_sfl = max(abs(par_sfl(:, (n_PS_kernel+2):end, :)), [], "all");
clim_all = max(clim_ori, clim_sfl);
clim_all = 2;

%% sorting
sorting_ranges = [1, A_T_border-0.5;A_T_border+0.5,T_P_border-0.5;T_P_border+0.5, N];
criterion = 1;
sort_idx = zeros(1, N);

criterion_mat = par_sig(:, (N*(criterion-1) + n_PS_kernel + 2):(N*criterion + n_PS_kernel + 1));
criterion_mat(isnan(criterion_mat))=0;
% get sorting index
for i=1:3
    sort_s = sorting_ranges(i, 1);
    sort_e = sorting_ranges(i, 2);

    if sorting == "count"
        s = criterion_mat~=0;
        inward = sum(s, 2).';
        outward = sum(s, 1);
        v = inward+outward;
        [~, sort_idx(sort_s:sort_e)] = sort(v(sort_s:sort_e),'descend');
    elseif sorting == "abs"
        s = abs(criterion_mat);
        inward = sum(s, 2).';
        outward = sum(s, 1);
        v = inward+outward;
        [~, sort_idx(sort_s:sort_e)] = sort(v(sort_s:sort_e), 'descend');
    elseif sorting == "square"
        s = criterion_mat.^2;
        inward = sum(s, 2).';
        outward = sum(s, 1);
        v = inward+outward;
        [~, sort_idx(sort_s:sort_e)] = sort(v(sort_s:sort_e), 'descend');
    elseif sorting == "max"
        s = abs(criterion_mat);
        inward = max(s, [], 2).';
        outward = max(s, [], 1);
        v = max(inward,outward);
        [~, sort_idx(sort_s:sort_e)] = sort(v(sort_s:sort_e), 'descend');
    elseif sorting == "idx"
        [~, sort_idx(sort_s:sort_e)] = sort(idx(sort_s:sort_e), 'ascend');
    end
    sort_idx(sort_s:sort_e) = sort_idx(sort_s:sort_e) + sort_s-1;
end


%% sorting
sort_s = sorting_ranges(i, 1);
sort_e = sorting_ranges(i, 2);
par_ori(:, :) = par_ori(sort_idx, :);
par_sig(:, :) = par_sig(sort_idx, :);
par_sfl(:, :, :) = par_sfl(sort_idx, :, :);
for j=1:n_conn_kernel
    shift = N*(j-1) + n_PS_kernel + 1;
    par_ori(:, shift+1:shift+N) = par_ori(:, shift+sort_idx);
    par_sig(:, shift+1:shift+N) = par_sig(:, shift+sort_idx);
    par_sfl(:, shift+1:shift+N, :) = par_sfl(:, shift+sort_idx, :);
end


%% save sorting index
sort_path = ['../GLM_data/', dataset_name, '/sortidx_', dataset_name, '_', ...
        int2str(session),'_', kernel_name, '.mat'];
save(sort_path, "sort_idx");


%% plotting
% kernel comparason
fig = figure("Visible", "off");
set(fig, 'PaperPosition', [0, 0, 12, 12]);
tiles = tiledlayout(3, 3);
cmap="jet";

kernel1 = par_ori(:, 4:3+N);
kernel2 = par_ori(:, 4+N:3+2*N);

% TODO: compare entries in kernel1 and kernel2.
% each subplot is an area.
% use A_T_border and T_P_border as area divide on both axis



fig_path = ['../figures/GLM/', dataset_name];
check_path(fig_path);
fig_file = [fig_path, '/GLMkernels_' dataset_name, '_', ...
        int2str(session), '_', kernel_name, '_', ...
        reg.name, '_', int2str(epoch), '_sorted.png'];
% exportgraphics(fig, fig_file);
print(fig, fig_file,'-dpng', '-r100');

end
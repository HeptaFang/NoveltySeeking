
% dataset_name = 'MuscimolPre_full'
dataset_name = 'SalinePost_cortex'
session_num = 5;
kernel_name = 'expDecay10';
reg.name = 'L2=5';
epoch = 2500;

for session = 1:session_num
    % load data
    data_path = ['../GLM_data/', dataset_name, '/GLMdata_', dataset_name, '_', ...
        int2str(session),'_', kernel_name, '_0.mat'];
    model_path = ['../GLM_model/', dataset_name, '/GLM_', dataset_name, '_', ...
        int2str(session), '_', kernel_name, '_0_', ...
        reg.name, '_', int2str(epoch)];
    border_path = ['../GLM_data/', dataset_name,'/borders_', dataset_name, '_', ...
        int2str(session),'.mat'];
    sort_path = ['../GLM_data/', dataset_name, '/sortidx_', dataset_name, '_', ...
        int2str(session),'_', kernel_name, '.mat'];
    raster_path = ['../GLM_data/', dataset_name, '/raster_', dataset_name, '_', ...
        int2str(session), '_0.mat'];

    load(data_path, 'raster');
    load(model_path, "model_par", "N");
    load(border_path, "borders");
    load(sort_path, "sort_idx");
    load(raster_path, "n_trial", "rasters", "spikes", "channel", "cell_area", "cell_id");

    % matrix and firing rate
    firing_rate = mean(raster, 2)*1000;
    J = model_par(:, 2:N+1);

    check_path(['../GLM_data/pack/', dataset_name])
    save_path = ['../GLM_data/pack/', dataset_name, '/connectivity_', dataset_name, '_', ...
        int2str(session), '_', kernel_name, '_0_', ...
        reg.name, '_', int2str(epoch)];
    save(save_path, "J", "firing_rate", "N", "spikes", "rasters", "n_trial",...
        "channel", "cell_area", "cell_id");

end

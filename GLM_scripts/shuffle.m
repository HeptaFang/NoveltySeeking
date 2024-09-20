function shuffle(dataset_name, session, shuffle_seed, across_trial)
%% shuffle spike trains and keep the firing rate.

% load original raster
raster_file = ['../GLM_data/', dataset_name, '/raster_', dataset_name, '_', int2str(session), '_0.mat'];
load(raster_file, "rasters", "n_trial", "trial_len");
firing_rates = load(raster_file).firing_rates;
n_raster = length(rasters);

rasters_shuffle = cell(1, n_raster);
rng(shuffle_seed);

if ~across_trial
    % shuffle within trial
    for i =1:n_raster
        raster = rasters{i};
        [N, B] = size(raster);
        raster_shuffle = zeros(N, B);
        for j = 1:N
            shuffle_idx = randperm(B);
            raster_shuffle(j, :) = raster(j, shuffle_idx);
        end
        rasters_shuffle{i} = raster_shuffle;
    end
else
    % shuffle across trial
    [N, ~] = size(rasters{1});
    full_raster = zeros(N, 0);
    for i=1:n_raster
        full_raster = [full_raster, rasters{i}];
    end
    [~, B] = size(full_raster);
    raster_shuffle = zeros(N, B);
    for j = 1:N
        shuffle_idx = randperm(B);
        raster_shuffle(j, :) = full_raster(j, shuffle_idx);
    end
    pointer=0;
    for i=1:n_raster
        Bi = size(rasters{i}, 2);
        rasters_shuffle{i} = raster_shuffle(:, pointer+1:pointer+Bi);
        pointer = pointer+Bi;
    end
end

raster_file_shuffle = ['../GLM_data/', dataset_name, '/raster_', ...
    dataset_name, '_', int2str(session),  '_', int2str(shuffle_seed), '.mat'];
rasters = rasters_shuffle;
save(raster_file_shuffle, "rasters", "firing_rates", "n_trial", "trial_len");
end
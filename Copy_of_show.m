pre = load('GLM_data/MuscimolPre_cortex/GLMdata_MuscimolPre_cortex_10_expDecay10_0.mat').raster;
post = load('GLM_data/MuscimolPost_cortex/GLMdata_MuscimolPost_cortex_10_expDecay10_0.mat').raster;

selected = cell(1, 0);

% selected = [selected, pre([14,51], 1:300000).'];
% selected = [selected, 0];
% selected = [selected, post([14,51], 1:300000).'];

selected = [selected, pre([8,21,23,9,10,11,25], :).'];
% selected = [selected, 0];
selected = [selected, post([8,21,23,9,10,11,25], :).'];


% align
maxlen = max(cellfun(@(arr) size(arr, 1), selected));
full = zeros(maxlen, 0);
for i=1:length(selected)
    padded = [selected{i};zeros(maxlen-size(selected{i}, 1), size(selected{i}, 2))];
    full = [full, padded];
end
full = full+(1:size(full, 2))-0.5;

% plot: L*N array
f = figure;
plot(full);
set(gca, 'YDir','reverse');
ylabel('cell id');
xlabel('time (ms)');
xline(0:3070:1135900);
pre = load('GLM_data/MuscimolPre_cortex/GLMdata_MuscimolPre_cortex_10_expDecay10_0.mat').raster;
post = load('GLM_data/MuscimolPost_cortex/GLMdata_MuscimolPost_cortex_10_expDecay10_0.mat').raster;

selected = cell(1, 0);

% selected = [selected, pre([14,51], 1:300000).'];
% selected = [selected, 0];
% selected = [selected, post([14,51], 1:300000).'];

selected = [selected, pre([8,21,23,9,10,11,25], :).'];
% selected = [selected, 0];
selected = [selected, post([8,21,23,9,10,11,25], :).'];


% align
maxlen = max(cellfun(@(arr) size(arr, 1), selected));
full = zeros(maxlen, 0);
for i=1:length(selected)
    padded = [selected{i};zeros(maxlen-size(selected{i}, 1), size(selected{i}, 2))];
    full = [full, padded];
end
full = full+(1:size(full, 2))-0.5;

% plot: L*N array
f = figure;
plot(full);
set(gca, 'YDir','reverse');
ylabel('cell id');
xlabel('time (ms)');
xline(0:3070:1135900);
pre = load('GLM_data/MuscimolPre_cortex/GLMdata_MuscimolPre_cortex_10_expDecay10_0.mat').raster;
post = load('GLM_data/MuscimolPost_cortex/GLMdata_MuscimolPost_cortex_10_expDecay10_0.mat').raster;

selected = cell(1, 0);

% selected = [selected, pre([14,51], 1:300000).'];
% selected = [selected, 0];
% selected = [selected, post([14,51], 1:300000).'];

selected = [selected, pre([8,21,23,9,10,11,25], :).'];
% selected = [selected, 0];
selected = [selected, post([8,21,23,9,10,11,25], :).'];


% align
maxlen = max(cellfun(@(arr) size(arr, 1), selected));
full = zeros(maxlen, 0);
for i=1:length(selected)
    padded = [selected{i};zeros(maxlen-size(selected{i}, 1), size(selected{i}, 2))];
    full = [full, padded];
end
full = full+(1:size(full, 2))-0.5;

% plot: L*N array
f = figure;
plot(full);
set(gca, 'YDir','reverse');
ylabel('cell id');
xlabel('time (ms)');
xline(0:3070:1135900);
pre = load('GLM_data/MuscimolPre_cortex/GLMdata_MuscimolPre_cortex_10_expDecay10_0.mat').raster;
post = load('GLM_data/MuscimolPost_cortex/GLMdata_MuscimolPost_cortex_10_expDecay10_0.mat').raster;

selected = cell(1, 0);

% selected = [selected, pre([14,51], 1:300000).'];
% selected = [selected, 0];
% selected = [selected, post([14,51], 1:300000).'];

selected = [selected, pre([8,21,23,9,10,11,25], :).'];
% selected = [selected, 0];
selected = [selected, post([8,21,23,9,10,11,25], :).'];


% align
maxlen = max(cellfun(@(arr) size(arr, 1), selected));
full = zeros(maxlen, 0);
for i=1:length(selected)
    padded = [selected{i};zeros(maxlen-size(selected{i}, 1), size(selected{i}, 2))];
    full = [full, padded];
end
full = full+(1:size(full, 2))-0.5;

% plot: L*N array
f = figure;
plot(full);
set(gca, 'YDir','reverse');
ylabel('cell id');
xlabel('time (ms)');
xline(0:3070:1135900);
pre = load('GLM_data/MuscimolPre_cortex/GLMdata_MuscimolPre_cortex_10_expDecay10_0.mat').raster;
post = load('GLM_data/MuscimolPost_cortex/GLMdata_MuscimolPost_cortex_10_expDecay10_0.mat').raster;

selected = cell(1, 0);

% selected = [selected, pre([14,51], 1:300000).'];
% selected = [selected, 0];
% selected = [selected, post([14,51], 1:300000).'];

selected = [selected, pre([8,21,23,9,10,11,25], :).'];
% selected = [selected, 0];
selected = [selected, post([8,21,23,9,10,11,25], :).'];


% align
maxlen = max(cellfun(@(arr) size(arr, 1), selected));
full = zeros(maxlen, 0);
for i=1:length(selected)
    padded = [selected{i};zeros(maxlen-size(selected{i}, 1), size(selected{i}, 2))];
    full = [full, padded];
end
full = full+(1:size(full, 2))-0.5;

% plot: L*N array
f = figure;
plot(full);
set(gca, 'YDir','reverse');
ylabel('cell id');
xlabel('time (ms)');
xline(0:3070:1135900);
pre = load('GLM_data/MuscimolPre_cortex/GLMdata_MuscimolPre_cortex_10_expDecay10_0.mat').raster;
post = load('GLM_data/MuscimolPost_cortex/GLMdata_MuscimolPost_cortex_10_expDecay10_0.mat').raster;

selected = cell(1, 0);

% selected = [selected, pre([14,51], 1:300000).'];
% selected = [selected, 0];
% selected = [selected, post([14,51], 1:300000).'];

selected = [selected, pre([8,21,23,9,10,11,25], :).'];
% selected = [selected, 0];
selected = [selected, post([8,21,23,9,10,11,25], :).'];


% align
maxlen = max(cellfun(@(arr) size(arr, 1), selected));
full = zeros(maxlen, 0);
for i=1:length(selected)
    padded = [selected{i};zeros(maxlen-size(selected{i}, 1), size(selected{i}, 2))];
    full = [full, padded];
end
full = full+(1:size(full, 2))-0.5;

% plot: L*N array
f = figure;
plot(full);
set(gca, 'YDir','reverse');
ylabel('cell id');
xlabel('time (ms)');
xline(0:3070:1135900);
pre = load('GLM_data/MuscimolPre_cortex/GLMdata_MuscimolPre_cortex_10_expDecay10_0.mat').raster;
post = load('GLM_data/MuscimolPost_cortex/GLMdata_MuscimolPost_cortex_10_expDecay10_0.mat').raster;

selected = cell(1, 0);

% selected = [selected, pre([14,51], 1:300000).'];
% selected = [selected, 0];
% selected = [selected, post([14,51], 1:300000).'];

selected = [selected, pre([8,21,23,9,10,11,25], :).'];
% selected = [selected, 0];
selected = [selected, post([8,21,23,9,10,11,25], :).'];


% align
maxlen = max(cellfun(@(arr) size(arr, 1), selected));
full = zeros(maxlen, 0);
for i=1:length(selected)
    padded = [selected{i};zeros(maxlen-size(selected{i}, 1), size(selected{i}, 2))];
    full = [full, padded];
end
full = full+(1:size(full, 2))-0.5;

% plot: L*N array
f = figure;
plot(full);
set(gca, 'YDir','reverse');
ylabel('cell id');
xlabel('time (ms)');
xline(0:3070:1135900);
pre = load('GLM_data/MuscimolPre_cortex/GLMdata_MuscimolPre_cortex_10_expDecay10_0.mat').raster;
post = load('GLM_data/MuscimolPost_cortex/GLMdata_MuscimolPost_cortex_10_expDecay10_0.mat').raster;

selected = cell(1, 0);

% selected = [selected, pre([14,51], 1:300000).'];
% selected = [selected, 0];
% selected = [selected, post([14,51], 1:300000).'];

selected = [selected, pre([8,21,23,9,10,11,25], :).'];
% selected = [selected, 0];
selected = [selected, post([8,21,23,9,10,11,25], :).'];


% align
maxlen = max(cellfun(@(arr) size(arr, 1), selected));
full = zeros(maxlen, 0);
for i=1:length(selected)
    padded = [selected{i};zeros(maxlen-size(selected{i}, 1), size(selected{i}, 2))];
    full = [full, padded];
end
full = full+(1:size(full, 2))-0.5;

% plot: L*N array
f = figure;
plot(full);
set(gca, 'YDir','reverse');
ylabel('cell id');
xlabel('time (ms)');
xline(0:3070:1135900);
pre = load('GLM_data/MuscimolPre_cortex/GLMdata_MuscimolPre_cortex_10_expDecay10_0.mat').raster;
post = load('GLM_data/MuscimolPost_cortex/GLMdata_MuscimolPost_cortex_10_expDecay10_0.mat').raster;

selected = cell(1, 0);

% selected = [selected, pre([14,51], 1:300000).'];
% selected = [selected, 0];
% selected = [selected, post([14,51], 1:300000).'];

selected = [selected, pre([8,21,23,9,10,11,25], :).'];
% selected = [selected, 0];
selected = [selected, post([8,21,23,9,10,11,25], :).'];


% align
maxlen = max(cellfun(@(arr) size(arr, 1), selected));
full = zeros(maxlen, 0);
for i=1:length(selected)
    padded = [selected{i};zeros(maxlen-size(selected{i}, 1), size(selected{i}, 2))];
    full = [full, padded];
end
full = full+(1:size(full, 2))-0.5;

% plot: L*N array
f = figure;
plot(full);
set(gca, 'YDir','reverse');
ylabel('cell id');
xlabel('time (ms)');
xline(0:3070:1135900);
pre = load('GLM_data/MuscimolPre_cortex/GLMdata_MuscimolPre_cortex_10_expDecay10_0.mat').raster;
post = load('GLM_data/MuscimolPost_cortex/GLMdata_MuscimolPost_cortex_10_expDecay10_0.mat').raster;

selected = cell(1, 0);

% selected = [selected, pre([14,51], 1:300000).'];
% selected = [selected, 0];
% selected = [selected, post([14,51], 1:300000).'];

selected = [selected, pre([8,21,23,9,10,11,25], :).'];
% selected = [selected, 0];
selected = [selected, post([8,21,23,9,10,11,25], :).'];


% align
maxlen = max(cellfun(@(arr) size(arr, 1), selected));
full = zeros(maxlen, 0);
for i=1:length(selected)
    padded = [selected{i};zeros(maxlen-size(selected{i}, 1), size(selected{i}, 2))];
    full = [full, padded];
end
full = full+(1:size(full, 2))-0.5;

% plot: L*N array
f = figure;
plot(full);
set(gca, 'YDir','reverse');
ylabel('cell id');
xlabel('time (ms)');
xline(0:3070:1135900);

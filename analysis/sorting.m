    
for session=1:10
    taskname = 'MuscimolPre_full';
    borderfile = ['../GLM_data/', taskname,'/borders_', taskname, '_', ...
    int2str(session),'.mat'];
    load(borderfile, "borders");
    channel_file = ['../GLM_data/', taskname, '/raster_', taskname, '_', int2str(session), ...
    '_0.mat'];
    load(channel_file, "channel", "N");
    borders = [1, round(borders+0.5, 0)];

    sort_idx = zeros(1, N);
    for i=1:3
        sort_s = borders(i);
        sort_e = borders(i+1)-1;
        [~, sort_idx(sort_s:sort_e)] = sort(channel(sort_s:sort_e), 'ascend');
        sort_idx(sort_s:sort_e) = sort_idx(sort_s:sort_e) + sort_s-1;
    end

    sort_file = ['../GLM_data/', taskname, '/sortIndex_', taskname, '_',...
        int2str(session), '.mat'];
    save(sort_file, "sort_idx");
end
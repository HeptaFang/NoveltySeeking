
%% generate poisson process
total=1000000;
on=10000;
off=total-on;

raster = [zeros(1, off), ones(1, on)];
raster = raster(randperm(total));

%% hold time: exponential distribution
all_isi=[];
count=0;
for i=1:total
    count=count+1;
    if raster(i)==1
        all_isi=[all_isi, count];
        count=0;
    end
end

figure;
histogram(all_isi, 0:1000,Normalization="probability",LineStyle="none");
hold on;

dist = (off/total).^(0:1000);
dist = dist/sum(dist);
plot(1:1001,dist);

%% fano factor, different time window
max_t = 10000;
fanos = zeros(1, max_t);
for t=1:max_t
    kernel = ones(1, t);
    counts = conv(raster, kernel, "valid");
    fanos(t) = var(counts)/mean(counts);
end
figure;
plot(fanos);


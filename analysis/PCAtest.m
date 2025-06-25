
clear; clc; close;

N = 64;
T = 100000;


myGauss = @(x, mu, sigma) exp(-((x-mu).^2)/(2*sigma^2));
% gaussian smoothing, sigma = 60ms, width = 300ms
kernel_x = -500:500;
kernel_y = myGauss(kernel_x, 0, 60);
smooth_kernel = kernel_y/sum(kernel_y); % normalize to 1

noise = randn(N, T)*0.25;
wave1 = sin(2*pi*0.01*(1:T));
wave2 = sin(4.3*pi*0.02*(1:T));

wave_complex1 = zeros(N, T);
wave_complex1(1:10, :) = repmat(wave1, 10, 1);
wave_complex1(11:20, :) = repmat(wave2, 10, 1);

activity1 = noise + wave_complex1*2;
activity1 = conv2(activity1, smooth_kernel, 'same');
activity1 = zscore(activity1, 0, 2);

activity2 = noise + wave2*1;
activity2 = conv2(activity2, smooth_kernel, 'same');
activity2 = zscore(activity2, 0, 2);

figure;
subplot(2,1,1);
plot(activity1(1:5, :)');

subplot(2,1,2);
plot(activity2(1:5, :)');


[coeff1,score1,latent1,tsquared1,explained1,mu1] = pca(activity1');
explained1_cumsum = cumsum(explained1);

[coeff2,score2,latent2,tsquared2,explained2,mu2] = pca(activity2');
explained2_cumsum = cumsum(explained2);

figure; hold on;
plot(explained1_cumsum);
plot(explained2_cumsum);



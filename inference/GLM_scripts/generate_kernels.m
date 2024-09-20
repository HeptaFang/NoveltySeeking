%% Util functions

% gaussian function
gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;

%% kernel generate 
% Each kernel set can contain multiple 'conn_kernel's (connection) 
% and 'PS_kernel's (post-spike). All kernels must have same length,
% add zeros to align all kernels.

% kernel file: "conn_kernels", "PS_kernels",
% "n_conn_kernel", "n_PS_kernel", "kernel_len"

% exponential decay kernel
tau1=10; % synaptic integration time constant (in ms)
T1=3*tau1; % cutoff on the sum for the neuron total input
tau_all1=0:T1;
tt_start=1+T1;
kernel = exp(-tau_all1./tau1);
kernel(1) = 0; % remove simultanuous corr
conn_kernels = {kernel};
PS_kernels = {};
n_conn_kernel=1;
n_PS_kernel=0;
kernel_len=T1+1;
save(['../GLM_data/', 'kernel_expDecay10.mat'], "conn_kernels", "PS_kernels", ...
    "n_conn_kernel", "n_PS_kernel", "kernel_len");

% gaussian kernel


% multi-kernel group



function GLM_multi_kernel(dataset_name, session, kernel_name, shuffle_seed, max_epoch, reg, log_level)
%% GLM inference.
%%%% required input : (from convolution)
%%%% dataset: "../GLM_data/[dataset_name]/
%%%%   GLMdata_[dataset_name]_[session]_[kernel_name][_shuffle].mat"
if nargin < 7
    log_level=2;
end
if nargin < 6
    reg.l1=0;
    reg.l2=0;
    reg.name='None';
end
if nargin < 5
    max_epoch=1000;
end


foldername=['../GLM_data/', dataset_name];
data_filename=[foldername, '/GLMdata_', dataset_name, '_', ...
    int2str(session), '_', kernel_name, '_', int2str(shuffle_seed), '.mat'];
load(data_filename, "raster", "predjs_conn", "predjs_PS", ...
    "conn_kernels", "PS_kernels", ...
    "n_conn_kernel", "n_PS_kernel", "kernel_len", "N", "B"); 

% Initial condition (can be better):
% h ~ (:, 1)
% P_ik ~ (:, k+1)
% J_ijk ~ (:, (N*(k-1) + n_PS_kernel + 2):(N*k + n_PS_kernel + 1))
par0=zeros(N, 1 + n_PS_kernel + N*n_conn_kernel); 

% Adam solver
beta1 = 0.9;
beta2 = 0.999;
lr = 1e-3 * sqrt(B/16); % large batch payoff
e = 1e-8;

m = zeros(size(par0));
v = zeros(size(par0));
par = par0;

% GPU acceleration
raster = gpuArray(raster);
predjs_PS = gpuArray(predjs_PS);
predjs_conn = gpuArray(predjs_conn);
par = gpuArray(par);
facts = factorial(raster);
logfacts = log(facts);
    
beta1_t = 1;
beta2_t = 1;
fprintf("Ready\n");
for epoch=1:max_epoch
    [minuslogL, grad_minuslogL, regLoss] = minuslogL_grad_fun(par,B,N, ...
        n_PS_kernel,n_conn_kernel,raster,predjs_PS,predjs_conn,logfacts,reg); 

    % Adam
    beta1_t = beta1_t * beta1;
    beta2_t = beta2_t * beta2;
    m = beta1*m + (1-beta1)*grad_minuslogL;
    v = beta2*v + (1-beta2)*grad_minuslogL.^2;
    mh = m/(1-beta1_t);
    vh = v/(1-beta2_t);
    
    % grad update
    par = par - lr*mh./(sqrt(vh+e));
    
    if log_level==2
        fprintf("Epoch %d/%d, -logL=%f, reg=%f, sparsity=%f\n", epoch, max_epoch, ...
            minuslogL, regLoss, sum(abs(par)>1e-2,"all")/numel(par));
    end
    
    % save model
    if mod(epoch, 100)==0
        if log_level==1
            fprintf("Epoch %d/%d, -logL=%f, reg=%f, sparsity=%f\n", epoch, max_epoch, ...
                minuslogL, regLoss, sum(abs(par)>1e-2,"all")/numel(par));
        end
        foldername = ['../GLM_model/', dataset_name];
        check_path(foldername);

        model_path = [foldername, '/GLM_', dataset_name, '_', ...
        int2str(session), '_', kernel_name, '_', int2str(shuffle_seed), '_', ...
        reg.name, '_', int2str(epoch)];

        [model_par, loss] = gather(par, minuslogL);
        save(model_path, 'model_par', 'loss', 'N', "reg", "kernel_len", ...
            "PS_kernels", "conn_kernels", "n_PS_kernel", "n_conn_kernel");
    end
end

  


end
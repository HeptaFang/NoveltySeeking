function p=GLM_test(par,B,N, ...
        n_PS_kernel,n_conn_kernel,raster,predjs_PS,predjs_conn,logfacts,reg)
% raster:(N, B), predj:(N, B, n_kernel)
% calculate p_it based on current parameters.
% return: p(N, B)

chunk_size_B = 200000;

hi=par(:, 1); % (N, 1)
Pik=par(:, 2:(n_PS_kernel+1)); % (N, n_PS)
Jijk=reshape(par(:, (n_PS_kernel+2):end), N, N, n_conn_kernel); % (N, N, n_conn)

htot_it = repmat(hi, 1, B); % (N, B)

%%%% tensor notations: 
%%%% <A>=elementwise product, 
%%%% [A]=contracted tensor product

% conn kernels
for k=1:n_conn_kernel
    htot_it = htot_it + Jijk(:, :, k) * predjs_conn(:, :, k); % (N, [N])*([N], B)=(N, B)
end
% post-spike kernels
for k=1:n_PS_kernel
    htot_it = htot_it + diag(Pik(:, k)) * predjs_PS(:, :, k); % (<N>)*(<N>, B)=(N, B)
end

% link function for poisson distribution
lambda_it=exp(htot_it);  % (N, B)
% minuslogL= -sum(raster.*log(lambda_it) - log(1+lambda_it) ,'all');

p = lambda_it./(1+lambda_it);

end


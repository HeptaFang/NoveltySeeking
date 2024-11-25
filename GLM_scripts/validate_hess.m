N = 10;
B = 1000;
tt_start = 1;
lambda_it = rand(B, N)+1;
predj = rand(N, B);

tic;
nb_par=N^2;
hess_minuslogL=zeros(nb_par,nb_par);
for ii=1:N
    for tt=tt_start:B
        lambda_prime_it=lambda_it(tt,ii)/(1+lambda_it(tt,ii));
        grad_htot_it=zeros(nb_par,1);
        grad_htot_it(ii,1)=1;
        predj_new=predj(:,tt);
        predj_new(ii)=[];    
        grad_htot_it(N+(ii-1)*(N-1)+(1:N-1),1)=predj_new;          
        for q1=1:nb_par
            for q2=1:nb_par
                hess_minuslogL(q1,q2)=hess_minuslogL(q1,q2)+lambda_prime_it*(1-lambda_prime_it)*grad_htot_it(q1,1)*grad_htot_it(q2,1);
            end
        end              
    end
end
% disp(hess_minuslogL);
err1 = sqrt(diag(inv(hess_minuslogL)));
toc;

tic;
% fprintf("calc hess\n");
n_PS_kernel = 0;
n_conn_kernel = 1;
predjs_PS = rand(N, B, n_PS_kernel);
predjs_conn = reshape(predj, N, B, 1);
err = zeros(N, 1 + n_PS_kernel + N*n_conn_kernel);
for i=1:N
    predj_PS_reshaped = permute(predjs_PS(i, :, :), [1, 3, 2]); % (1, n_PS, B)
    predj_PS_reshaped = reshape(predj_PS_reshaped, 1*n_PS_kernel, B); % (n_PS, B)
    predj_conn_reshaped = permute(predjs_conn, [1, 3, 2]); % (N, n_conn, B)
    predj_conn_reshaped = reshape(predj_conn_reshaped, N*n_conn_kernel, B); % (N*n_conn, B)
    % remove self-connections
    for k=1:n_conn_kernel
        predj_conn_reshaped(i + (k-1)*(N-1), :) = [];
    end

    grad_htot_it = [ones(1, B); predj_PS_reshaped; predj_conn_reshaped]; % (1 + n_PS + (N-1)*n_conn, B)
    lambda_factor_it = lambda_it(:, i)./((1+lambda_it(:, i)).^2); % (B, 1)
    lambda_factor_it_reshaped = repmat(lambda_factor_it.', 1 + n_PS_kernel + (N-1)*n_conn_kernel, 1); % (1 + n_PS + (N-1)*n_conn, B)
    hess = tensorprod(grad_htot_it, grad_htot_it.*lambda_factor_it_reshaped, 2, 2); % (1 + n_PS + (N-1)*n_conn, 1 + n_PS + (N-1)*n_conn)
    err_i = sqrt(diag(inv(hess))); % (1 + n_PS + (N-1)*n_conn)

    for k = 1:n_conn_kernel
        % put back zero to self-connections
        err_i = [err_i(1:(1+n_PS_kernel+(k-1)*N+i-1)); 0; err_i((1+n_PS_kernel+(k-1)*N+i):end)]; % (1 + n_PS + N*n_conn)
    end
    err(i, :) = err_i;
end
toc;
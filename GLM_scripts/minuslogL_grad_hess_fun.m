function [loss, grad, err]=minuslogL_grad_hess_fun(par,B,N, ...
        n_PS_kernel,n_conn_kernel,raster,predjs_PS,predjs_conn,logfacts,reg)
% raster:(N, B), predj:(N, B, n_kernel)

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

L_it = -raster.*log(lambda_it)+log(1+lambda_it);            % -log(L_it) for Bernoulli
% L_it = lambda_it - (raster.*log(lambda_it)-logfacts);       % -log(L_it) for Poisson

minuslogL = sum(L_it, "all");
loss.minuslogL = minuslogL;

RegLoss = 0;
if reg.l1>0 
    RegLoss = RegLoss + reg.l1*sum(abs(par(:, 2:end)), "all");
end
if reg.l2>0
    RegLoss = RegLoss + reg.l2*sum(par(:, 2:end).^2, "all");
end
loss.reg = RegLoss;
loss.total = loss.minuslogL + RegLoss;

if nargout > 1     
    % fprintf("calc grad\n");
    % todo: fix the grad calc
    % grad_minuslogL=zeros(N, 1 + n_PS_kernel + N*n_conn_kernel);
    
    dLdh_it = lambda_it./(1+lambda_it)-raster;  % d(minuslogL)/d(htot_it) for Bernoulli, (N, B)
    % dLdh_it = lambda_it-raster;                 % d(minuslogL)/d(htot_it) for Poisson, (N, B)

    dLdh_i = sum(dLdh_it, 2);% (N, [B])=(N, 1)
    dLdP_ik = sum(dLdh_it.*predjs_PS(:, :, :), 2);% (<N>, [B])*(<N>, [B], n_PS)=(N, n_PS)
    dLdJ_ijk = tensorprod(dLdh_it, predjs_conn, 2, 2); % (N, [B])*(N, [B], n_conn)=(N, N, n_conn)

    predj_conn_reshaped = permute(predjs_conn, [1, 3, 2]); % (N, n_conn, B)
    predj_conn_reshaped = reshape(predj_conn_reshaped, 1, N*n_conn_kernel, B); % (1, N*n_conn, B)
    predj_conn_reshaped = repmat(predj_conn_reshaped, N, 1, 1); % (N, N*n_conn, B)
    grad_htot_it = [ones(N, B, 1), predjs_PS, predj_conn_reshaped]; % (N, 1 + n_PS + N*n_conn, B)

    grad_minuslogL = sum(reshape(dLdh_it, N, 1, B).*grad_htot_it, 3); % (<N>, [B])*(<N>, 1 + n_PS + N*n_conn, [B])=(N, 1 + n_PS + N*n_conn)
    % eliminate self-connections in conn kernels
    for i=1:N
        for k=1:n_conn_kernel
            grad_minuslogL(i, i + n_PS_kernel + (k-1)*N + 1)=0;
        end
    end
    grad.minuslogL = grad_minuslogL;

    % grad_minuslogL(:, 1) = dLdh_i;
    % grad_minuslogL(:, 2:(n_PS_kernel+1)) = dLdP_ik;
    % grad_minuslogL(:, (n_PS_kernel+2):end) = reshape(dLdJ_ijk, N, N*n_conn_kernel);
    
    % regularizations
    % if reg.l1>0 
    %     grad_minuslogL = grad_minuslogL + [zeros(N, 1), reg.l1*sign(par(:, 2:end))];
    % end
    % if reg.l2>0
    %     grad_minuslogL = grad_minuslogL + [zeros(N, 1), reg.l2*2*par(:, 2:end)];
    % end
    % if reg.l1>0 
    %     grad_minuslogL = grad_minuslogL + [zeros(N, 1), (1/(reg.l1*B))*sign(par(:, 2:end))];
    % end
    % if reg.l2>0
    %     grad_minuslogL = grad_minuslogL + [zeros(N, 1), (1/(2*reg.l2*B))*2*par(:, 2:end)];
    % end
    grad_reg = zeros(N, 1 + n_PS_kernel + N*n_conn_kernel);
    if reg.l1>0 
        grad_reg = grad_reg + [zeros(N, 1), reg.l1*sign(par(:, 2:end))];
    end
    if reg.l2>0
        grad_reg = grad_reg + [zeros(N, 1), reg.l2*2*par(:, 2:end)];
    end
    
    % eliminate self-connections in conn kernels
    for i=1:N
        for k=1:n_conn_kernel
            grad_reg(i, i + n_PS_kernel + (k-1)*N + 1)=0;
        end
    end

    grad.reg = grad_reg;
    grad.total = grad_minuslogL + grad_reg;

end

% Hessian matrix error (N, 1 + n_PS_kernel + N*n_conn_kernel)
if nargout > 2 
    % fprintf("calc hess\n");
    err = zeros(N, 1 + n_PS_kernel + N*n_conn_kernel);
    for i=1:N
        hess = tensorprod(grad_htot_it(i, :, :), grad_htot_it(i, :, :), [1, 3], [1, 3]); % (1 + n_PS + N*n_conn, 1 + n_PS + N*n_conn)
        err(i, :) = sqrt(diag(inv(hess)));
    end
end


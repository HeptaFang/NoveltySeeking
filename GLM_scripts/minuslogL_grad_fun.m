function [minuslogL,grad_minuslogL,RegLoss]=minuslogL_grad_fun(par,B,N, ...
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

if nargout > 1     
    % fprintf("calc grad\n");
    % todo: fix the grad calc
    grad_minuslogL=zeros(N, 1 + n_PS_kernel + N*n_conn_kernel);
    
    dLdh_it = lambda_it./(1+lambda_it)-raster;  % d(minuslogL)/d(htot_it) for Bernoulli
    % dLdh_it = lambda_it-raster;                 % d(minuslogL)/d(htot_it) for Poisson

    grad_hi = sum(dLdh_it, 2);% (N, [B])=(N, 1)
    grad_Pik = sum(dLdh_it.*predjs_PS(:, :, :), 2);% (<N>, [B])*(<N>, [B], n_PS)=(N, n_PS)
    grad_Jijk = tensorprod(dLdh_it, predjs_conn, 2, 2); % (N, [B])*(N, [B], n_conn)=(N, N, n_conn)

    grad_minuslogL(:, 1) = grad_hi;
    grad_minuslogL(:, 2:(n_PS_kernel+1)) = grad_Pik;
    grad_minuslogL(:, (n_PS_kernel+2):end) = reshape(grad_Jijk, N, N*n_conn_kernel);
    
    % regularizations
    if reg.l1>0 
        grad_minuslogL = grad_minuslogL + [zeros(N, 1), reg.l1*sign(par(:, 2:end))];
    end
    if reg.l2>0
        grad_minuslogL = grad_minuslogL + [zeros(N, 1), reg.l2*2*par(:, 2:end)];
    end

    % eliminate self-connections in conn kernels
    for i=1:N
        for k=1:n_conn_kernel
            grad_minuslogL(i, i + n_PS_kernel + (k-1)*N + 1)=0;
        end
    end

end


if nargout > 2 
    RegLoss = 0;
    if reg.l1>0 
        RegLoss = RegLoss + reg.l1*sum(abs(par(:, 2:end)), "all");
    end
    if reg.l2>0
        RegLoss = RegLoss + reg.l2*sum(par(:, 2:end).^2, "all");
    end
end



end


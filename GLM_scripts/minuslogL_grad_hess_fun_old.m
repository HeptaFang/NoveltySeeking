
function [minuslogL,grad_minuslogL,hess_minuslogL]=minuslogL_grad_hess_fun(par_vector,B,N,tt_start,raster,predj)

lambda_it=zeros(B,N);
for ii=1:N
    hi=par_vector(ii,1);
    Jij=par_vector(N+(ii-1)*(N-1)+[1:N-1],1)';  
    for tt=tt_start:B
        predj_new=predj(:,tt);
        predj_new(ii)=[];         
        htot_it=hi + Jij*predj_new;
        lambda_it(tt,ii)=exp(htot_it);
    end  
end
minuslogL= sum( - raster([tt_start:B],:).*log(lambda_it([tt_start:B],:)) + log(1+lambda_it([tt_start:B],:)) ,'all');


if nargout > 1     
    nb_par=N^2;
    grad_minuslogL=zeros(nb_par,1);
    for ii=1:N
        for tt=tt_start:B
            grad_htot_it=zeros(nb_par,1);
            grad_htot_it(ii,1)=1;
            predj_new=predj(:,tt);
            predj_new(ii)=[];    
            grad_htot_it(N+(ii-1)*(N-1)+[1:N-1],1)=predj_new;
            grad_minuslogL=grad_minuslogL-(raster(tt,ii)-(lambda_it(tt,ii)/(1+lambda_it(tt,ii))))*grad_htot_it;
        end
    end
end


if nargout > 2   
    nb_par=N^2;
    hess_minuslogL=zeros(nb_par,nb_par);
    for ii=1:N
        for tt=tt_start:B
            lambda_prime_it=lambda_it(tt,ii)/(1+lambda_it(tt,ii));
            grad_htot_it=zeros(nb_par,1);
            grad_htot_it(ii,1)=1;
            predj_new=predj(:,tt);
            predj_new(ii)=[];    
            grad_htot_it(N+(ii-1)*(N-1)+[1:N-1],1)=predj_new;          
            for q1=1:nb_par
                for q2=1:nb_par
                    hess_minuslogL(q1,q2)=hess_minuslogL(q1,q2)+lambda_prime_it*(1-lambda_prime_it)*grad_htot_it(q1,1)*grad_htot_it(q2,1);
                end
            end              
        end
    end
end



end


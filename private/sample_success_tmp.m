clear;
clc;
close;
m =100; n=100; r=3;  sample_budget = 2000; Alpha = [0.3]; Beta = [0.05]; 
Gama = [0.85];
num_iter = 50;  CVX_or_SVT = 0;
result_ = []; gamaInd = 1; alphaInd = 1; betaInd = 1;
for beta = Beta
    alphaInd = 1;
    for alpha = Alpha
        gamaInd = 1;
        for gama = Gama
            result = zeros(1,num_iter);
            for iter = 1:num_iter
                X0 = randn(m,r)*randn(r,n);
                if alpha ~=0
                    for pg=1:1
                        [u, s, v] = svd(X0);
                        u = u(:,1:r);
                        v = v(:,1:r);
                        s = s(1:r,1:r);
                        pgf1 = randi([0,m-1],1,1);
                        D1 = zeros(m,n);
                        for i =1:m
                            D1(i,i) = (abs(i-pgf1)+1)^(-alpha);
                        end
                        pgf2 = randi([0,m-1],1,1);
                        D2 = zeros(m,n);
                        for i =1:m
                            D2(i,i) = (abs(i-pgf2)+1)^(-alpha);
                        end
                        X0 = (D1 *u*s*v'* D2);
                        X0 = X0/norm(X0,'fro');
                    end
                end
                
                X0 = X0 / norm(X0,'fro');
                [U0, S0, V0] = svd(X0);
                num_first = floor(beta * sample_budget);
                random_numbers = randperm(m*n);
                omega = random_numbers(1:num_first);
                omega_diff = omega;
                if CVX_or_SVT==0
                    A = zeros(n,m);
                    A(omega) = 1;
                    y= vec(A.*X0);
                    cvx_begin quiet
                    variable Xt(m,n);
                    minimize norm_nuc(Xt);
                    subject to
                    y == vec(A.*Xt);
                    cvx_end
                else
                    [Ut,St,Vt,~] = Test_SVT(m,n,r,length(omega),omega,X0);
                    Xt = Ut*St*Vt';
                end
                error = norm( Xt - X0,'fro') / norm(X0,'fro');
                error_diff = 1; error_tereshold = 0.01;
                
                if error_diff > error_tereshold
                    [Xt_diff, omega_diff] = fun_diff_X0_Xn(X0, omega_diff, r, sample_budget - num_first, ceil((1-gama)*num_first),1);
                    error_diff = norm( Xt_diff - X0,'fro') / norm(X0,'fro');
                end
                result(iter) = error_diff;
                
                X = sprintf('error: %f  err_diff: %f  num: %d  gama: %f  alpha: %f  beta: %f  iter: %d'...
                    ,error,error_diff,length(omega_diff),gama,alpha,beta,iter);
                disp(X);
            end
            mean(result')
            result_(betaInd,alphaInd,gamaInd) = mean(result);
            X = sprintf('Beta_plot m_n_%d num_sample_ %d',m,sample_budget);
            save(char(X));
            gamaInd = gamaInd + 1;
        end
        alphaInd = alphaInd + 1;
    end
    betaInd = betaInd + 1;
end
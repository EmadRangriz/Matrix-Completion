function [Xt, omega] = fun_diff_X0_Xn(X0, omega,r, num_new_ind, diff,take_max)
m = length(X0(:,1));
n = length(X0(1,:));
omega = omega(randperm(length(omega)));
if diff ~= 0 && diff ~= length(omega)
    
        A = zeros(n,m);
        A(omega(1:end-diff)) = 1;
        y= vec(A.*X0);
        cvx_begin quiet
        variable XBn(m,n);
        minimize norm_nuc(XBn);
        subject to
        y == vec(A.*XBn);
        cvx_end
        
        A = zeros(n,m);
        A(omega) = 1;
        y= vec(A.*X0);
        cvx_begin quiet
        variable Xn(m,n);
        minimize norm_nuc(Xn);
        subject to
        y == vec(A.*Xn);
        cvx_end
        Xd = abs(Xn-XBn);
elseif diff==0
    Xd = ones(m,n);
elseif diff==length(omega)
    A = zeros(n,m);
    A(omega) = 1;
    y= vec(A.*X0);
    cvx_begin quiet
    variable Xd(m,n);
    minimize norm_nuc(Xd);
    subject to
    y == vec(A.*Xd);
    cvx_end
    Xd = abs(Xd);
end

cumsumed = cumsum(vec(Xd));
cumsumed = cumsumed/max(cumsumed);
j=0;
index = [];
if take_max ==1
    chInd = randperm(m*n);
    Xd(1:m*n)=Xd(chInd);
    [soF,indF] = sort(vec(Xd),'descend');
    for i = indF'
        if isempty(find(i == index, 1)) && isempty(find(i == omega, 1))
            index(j+1) = i;
            j = j+1;
        end
        if j >= num_new_ind
            break;
        end
    end
    index = chInd(index);
    Xd(chInd)=Xd(1:m*n);
else
    while length(index) < num_new_ind
        random_num = rand(1);
        for i = 1:m*n
            if random_num <= cumsumed(i)
                if isempty(find(i == index, 1)) && isempty(find(i == omega, 1))
                    index(j+1) = i;
                    j = j+1;
                end
                break;
            end
        end
    end
end
omega = [omega,index];

A = zeros(n,m);
A(omega) = 1;
y= vec(A.*X0);
cvx_begin quiet
variable Xt(m,n);
minimize norm_nuc(Xt);
subject to
y == vec(A.*Xt);
cvx_end

end
function [w] = clls(x1, x2, alpha)
% CLLS 

y = [x1, x2; ones(1, size(x1, 2) + size(x2, 2))];
N = size(y, 2);
L = size(y, 1);
c = [ones(1, size(x1, 2)), -1*ones(1, size(x2, 2))];

termo1 = zeros(L,L);
termo2 = zeros(L,1);
for k = 1:N
    termo1 = termo1 + y(:,k)*y(:,k)' + alpha*eye(L);
    termo2 = termo2 + c(k)*y(:,k);
end

w = (inv(termo1)*termo2)';

return
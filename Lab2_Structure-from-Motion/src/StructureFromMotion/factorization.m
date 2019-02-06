function [M, S] = factorization(D, eliminate_affinity)

assert(~any(any(isnan(D))));

[m, n] = size(D);
assert (mod(m, 2) == 0);
m = m / 2;

center = mean(D, 2);

D_n = D - center;

[U, W, V] = svd(D_n);

W = W(1:3, 1:3);
U = U(:, 1:3);
V = V(:, 1:3);

W_sqrt = sqrt(W);

M = U * W_sqrt;
S = W_sqrt * V';
  
% M = U;
% S = W*V';



%affinity ambiguity elimination
if eliminate_affinity
    A = zeros(4*m, 9);
    I = zeros(2*m, 2);
    
    
    %See "Compatibility with Kronecker products" in https://en.wikipedia.org/wiki/Vectorization_(mathematics)
    for i=1:m
        A_i = M(2*i-1:2*i, :);
        p = kron(A_i, A_i);
        A(4*i-3:4*i, :) = p;
        I(2*i-1:2*i, :) = eye(2);
    end

    I = reshape(I', [], 1);

    L_vec = I \ A;
    
    L = reshape(L_vec, 3, 3);
    
    C = chol(L);

    M = M*C;
    S = inv(C)*S;
end

% figure(11);
% scatter3(S(1, :)', S(2, :)', S(3, :)', 4, 'filled');

end


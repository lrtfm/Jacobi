function [ x, omega ] = JacobiGauss( alpha, beta, N )
    %JACOBIGAUSS Summary of this function goes here
    %   Detailed explanation goes here
    j = 0:N;
    T = 2*j + alpha + beta;
    D0 = (beta^2 - alpha^2)./(T .* (T + 2));
    j = 1:N;
    T = 2*j + alpha + beta;
    D1 = 4 * j .* (j + alpha) .* (j + beta) .* ( j + alpha + beta) ./...
        ((T - 1) .* (T.^2) .*(T+1));
    D1 = D1.^(1/2);
    A = diag(D0) + diag(D1, 1) + diag(D1, -1);
    if alpha == 0 && beta == 0
        A(1, 1) = 0;
    end
    [V, D] = eig(A);   % A 是是对成矩阵时， 这里的 V 是 orthonormal 的
    r = 2^(alpha + beta + 1)*gamma(alpha + 1)* gamma(beta + 1)/gamma(alpha + beta + 2);
    x = diag(D, 0);
    omega = (V(1, :).^2 * r)';
end


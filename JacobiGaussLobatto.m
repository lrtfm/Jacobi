function [ x, omega ] = JacobiGaussLobatto( alpha, beta, N )
    % JacobiGaussLobatto
    
    [ y, w ] = Jacobi.JacobiGauss( alpha+1, beta+1, N - 2 );
    x = [ -1; y; 1];
    g0 = 2*gammaln(beta + 1) + gammaln(N) + gammaln(N + alpha + 1) -...
        gammaln(N + beta + 1) - gammaln(N + alpha + beta + 2);
    omega0 = 2^(alpha + beta + 1) * (beta + 1) * exp(g0);
    gN = 2*gammaln(alpha + 1) + gammaln(N) + gammaln(N + beta + 1) -...
        gammaln(N + alpha + 1) - gammaln(N + alpha + beta + 2);
    omegaN = 2^(alpha + beta + 1) * (alpha + 1) * exp(gN);
    omega = [omega0; w ./ (1 - y.^2); omegaN];
end
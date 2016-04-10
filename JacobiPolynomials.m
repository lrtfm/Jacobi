function J = JacobiPolynomials(alpha, beta, N, x)
%     J = 0;
%     for k = 0:N;
%         J = J + 1/(factorial(k)*factorial(N-k)) *...
%             gamma(N + k + alpha + beta + 1) / gamma(k + alpha + 1) * ((x - 1)/2).^k;
%     end
%     J = gamma(N + alpha + 1) / gamma(N + alpha + beta + 1) * J;
    if N == 0
        J = ones(size(x));
        return
    elseif N == 1
        J = (alpha + beta + 2)/2*x + (alpha - beta)/2;
        return
    end
    
    PJ = ones(size(x)); J = (alpha + beta + 2)/2*x + (alpha - beta)/2;
    for n = 1:(N-1)
        PPJ = PJ; 
        PJ = J;
        an = (2*n + alpha + beta + 1) * (2*n + alpha + beta + 2) ...
            / (2 * (n + 1) * (n + alpha + beta + 1));
        bn = (beta^2 - alpha^2)*(2*n + alpha + beta + 1) ...
            / ( 2*(n + 1) * (n + alpha + beta + 1) * (2*n + alpha + beta));
        cn = (n + alpha)*(n + beta)*(2*n + alpha + beta + 2) / ( (n + 1) ...
            * (n + alpha + beta + 1) * (2*n + alpha + beta));
        J = (an*x - bn) .* PJ - cn * PPJ;
    end
end
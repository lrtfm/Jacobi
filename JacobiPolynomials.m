function [JN, J] = JacobiPolynomials(alpha, beta, N, x)
    x = x(:);
    J = zeros(length(x), N+1);
    for i = 0:N
        if i == 0
            J(:, i+1) = ones(length(x),1);
        elseif i == 1
            J(:, i+1) = (alpha + beta + 2)/2*x + (alpha - beta)/2;
        else
            n = i - 1;
            an = (2*n + alpha + beta + 1) * (2*n + alpha + beta + 2) ...
                / (2 * (n + 1) * (n + alpha + beta + 1));
            bn = (beta^2 - alpha^2)*(2*n + alpha + beta + 1) ...
                / ( 2*(n + 1) * (n + alpha + beta + 1) * (2*n + alpha + beta));
            cn = (n + alpha)*(n + beta)*(2*n + alpha + beta + 2) / ( (n + 1) ...
                * (n + alpha + beta + 1) * (2*n + alpha + beta));
            J(:, i+1) = (an*x - bn) .* J(:, i) - cn * J(:, i - 1);
        end
    end
    JN = J(:, end);
end

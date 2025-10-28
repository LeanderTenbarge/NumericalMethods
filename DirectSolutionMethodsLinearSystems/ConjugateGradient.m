% Leander Tenbarge
% Conjugate Gradient (Kylov subspace method)


function [u] = ConjugateGradient(A,c,tol,maxIter)
    if A' ~= A 
        error("The Matricies is Non-Symmemtric")
    end
    [~, p] = chol(A);
    if p ~= 0
        error('Matrix A must be positive-definite');
    end
    n = size(c,1);
    u = zeros(n,1);
    r = c;
    d = c; 
    i = 1;
    while max(r) > tol && i < maxIter
        rTemp = r;
        dTemp = d;
        alpha = ((rTemp')*(rTemp))/((dTemp') * (A) * (dTemp));
        u = u + alpha * dTemp;
        r = rTemp - alpha * (A) * dTemp;
        beta = ((r') * (r))/((rTemp') * (rTemp));
        d = r + beta * dTemp;
        i = i + 1;
    end    
end
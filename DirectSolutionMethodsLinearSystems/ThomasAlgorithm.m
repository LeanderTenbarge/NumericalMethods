% Leander Tenbarge
% Thomas Algorithm, Anderson CFD, page 151.

function [U] = ThomasAlgorithm(A,b)
    NJ = size(A,1);
    U = zeros(NJ,1);
    % Putting it into upper triangular form:
    for i = 2:NJ
        A(i,i) = A(i,i) - ((A(i,i-1))/(A(i-1,i-1))) * A(i-1,i);
        b(i) = b(i) - ((A(i,i-1))/(A(i-1,i-1))) * (b(i-1));
    end
    
    % Computing the unknowns from back substitution:
    U(NJ) = b(NJ)/A(NJ,NJ);
    for i = 2:NJ
        k = ( NJ + 1 - i);
        U(k) = (b(k) - A(k,k+1) * U(k+1))/(A(k,k));
    end
end
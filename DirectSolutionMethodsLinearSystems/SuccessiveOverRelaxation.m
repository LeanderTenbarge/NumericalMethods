% Leander Tenbarge
% Successive Over Relaxation method, Anderson CFD, page 152:

function [U] = SuccessiveOverRelaxation(A,b,tol,omega)

    % Processing the input:
    N = size(A,1);
    U = zeros(N,1);
    Residuals = ones(N,1);
    GreaterThan = false;
    
    % Checking the Stability Factor:
    for i = 1:N
        Magnitude = A(i,i);
        Sum = 0;
        for j = 1:N
            if j ~= i
                Sum = Sum + abs(A(i,j));
            end
        end

        % Checking the Stability Criterion
        if Magnitude < Sum
           error("The Summation of the Magnitude of the Off-Diagonals is larger than the Diagonal")
        end

        if Magnitude > Sum
            GreaterThan = true;
        end
    end

    if ~GreaterThan
        error("At no point the Summation of the Magnitude of the Diagonals is larger than the Off-Diagonal")
    end
    

    % Solving the System:
    while max(abs(Residuals)) > tol
        
        % Iteration loop:
        for i = 1:N
            % Initializing/Refreshing the solution variables:
            secondary = 0;
            tertiary = 0;

            % Primary Loop
            primary = b(i);
            
            % Secondary Loop:
            for j = 1:i-1
                secondary = secondary + A(i,j)*U(j);
            end

            % Tertiary Loop:
            for j = 1+i:N
                tertiary = tertiary + A(i,j)*U(j);
            end

            % Determining the solution:
            Ut = U(i);
            U(i) = (1-omega) * Ut+ (omega/A(i,i))* (primary - secondary - tertiary);
            Residuals(i) = (U(i) - Ut);
        end
    end
end


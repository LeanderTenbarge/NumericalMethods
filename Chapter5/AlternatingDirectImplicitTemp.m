% Leander Tenbarge: Hoffman CFD Chapter 5, Numerical Methods for Elliptic partial differential equations:
% Alternating Direct Implicit, Temperature Distribution:
% Page 175 - 176:

clear all

% System Setup as listed in page 155:
width = 1; 
height = 1; 
numPy = 150; % Number points in the width:
numPx = 150; % Number of points in the height:
dx = width/numPy;
dy = height/numPx;
n = 5000;
beta = dx/dy;
Tupper = 0;   % rankine
Tleft  = 0;   % rankine
Tright = 100; % rankine
Tlower = 100; % rankine
alpha = -2*(1+ beta^2); 
omega = 1;

% Implementing boundary conditions:
T = zeros(numPx,numPy); 
T(:,1) = Tleft;
T(1,:) = Tupper;
T(:,end) = Tright;
T(end,:) = Tlower;

Thalf = zeros(numPx,numPy); 
Thalf(:,1) = Tleft;
Thalf(1,:) = Tupper;
Thalf(:,end) = Tright;
Thalf(end,:) = Tlower;

% Setting up the X sweep:
A1 = spdiags([1,alpha,1],[-1,0,1],numPx-2,numPx-2);

% Setting up th Y sweep:
A2 = spdiags([beta^2,alpha,beta^2],[-1,0,1],numPy-2,numPy-2);

for iter = 1:n

    % Evaluating the First equation:
    for j = 2:numPy-1

        % Enforcing the Boundary Conditions:
        b1 = zeros(numPy-2,1);
        b1(1) = -T(1,j);
        b1(end) = -T(end,j);
        
        for i = 1:numPx - 2
            b1(i) = b1(i) + (-beta^2) * (T(i+1,j+1) + T(i+1,j-1));
        end

        % Solving the Tridiagional System of equations:
        Thalf(2:numPx-1,j) = A1\b1;

    end

    % Evaluating the Second equation:
    for i = 2:numPx-1

        % Enforcing the Boundary Conditions:
         b2 = zeros(numPy-2,1);
         b2(1) = -(beta^2)*Thalf(i,1);
         b2(end) = -(beta^2)*Thalf(i,end);

        for j = 1:numPy - 2
            b2(j) = b2(j) - Thalf(i+1,j+1) - Thalf(i-1,j+1);
        end
        
        % Solving the tridiagonal system of equations
        T(i,2:numPy-1) = A2\b2;
    end 
end

% Plotting the results
contour(T);        
colorbar;             
colormap('jet');        
title('ADI Formulation, Steady State Heat equation');
xlabel('X-axis');
ylabel('Y-axis');
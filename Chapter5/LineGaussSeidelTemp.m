% Leander Tenbarge, Hoffman CFD Chapter 5, Numerical Methods for Elliptic partial differential equations:
% Line Gauss Seidel Method:
% Page 153-155:
clear all
% System Setup as listed in page 155:
width = 1; 
height = 1; 
numPW = 50; % Number points in the width:
numPH = 50; % Number of points in the height:
dx = width/numPW;
dy = height/numPH;
n = 500;
beta = dx/dy;
Tupper = 0;   % rankine
Tleft  = 0;   % rankine
Tright = 100; % rankine
Tlower = 100; % rankine
alpha = -2*(1+ beta^2);

% Implementing boundary conditions:
T = zeros(numPH,numPW); 
T(:,1) = Tleft;
T(1,:) = Tupper;
T(:,end) = Tright;
T(end,:) = Tlower;

% Implementing the Line Gauss - Seidel method:


% Creating the A matricies:
A = spdiags([1,alpha,1],[-1,0,1],numPH-2,numPH-2);

% Looping over time:
for time = 1:n
    for j = 2:numPW-1

        % Generating the Matricies:
        B = zeros(numPH-2,1);

        % Enforcing the Boundary conditions:
        B(1) = -T(1,j);
        B(numPH-2) = -T(end,j);
        
        % Creating the rest of the left hand side:
        for i = 1:numPH-2
            B(i) = B(i) + (-beta^2)*(T(i+1,j+1) + T(i+1,j-1));
        end
    end
end

% Plotting the results
imagesc(T);        
colorbar;             
colormap('jet');        
title('Line Gauss Siedel, Steady State Heat equation');
xlabel('X-axis');
ylabel('Y-axis');

% Leander Tenbarge:
% Point Gauss - Seidel for Temperature Distribution:
% Page 175 - 176:
clear all
% System Setup as listed in page 155:
width = 1; 
height = 1; 
numPW = 500; % Number points in the width:
numPH = 500; % Number of points in the height:
dx = width/numPW;
dy = height/numPH;
n = 100;
beta = dx/dy;
Tupper = 0;   % rankine
Tleft  = 0;   % rankine
Tright = 100; % rankine
Tlower = 100; % rankine

% Implementing boundary conditions:
T = zeros(numPH,numPW); 
T(:,1) = Tleft;
T(1,:) = Tupper;
T(:,end) = Tright;
T(end,:) = Tlower;


% Running over the Gauss - Seidel algorithm:

for time = 1:n

    % Boundary Conditions for temp:
    Ttemp = zeros(numPH,numPW);
    Ttemp(:,1) = Tleft;
    Ttemp(1,:) = Tupper;
    Ttemp(:,end) = Tright;
    Ttemp(end,:) = Tlower;

    for j = 2:numPW-1
        for i = 2:numPH-1
            Ttemp(i,j) = (1/(2*(1+beta^2))) * (T(i+1,j) + Ttemp(i-1,j) + (beta ^ 2)*(T(i,j+1) + Ttemp(i,j-1)));
        end
    end
    
    % Implementing the solution to the main matricies:
    T = Ttemp;
end


% Plotting the results
imagesc(T);        
colorbar;             
colormap('jet');        
title('Gauss Siedel, Steady State Heat equation');
xlabel('X-axis');
ylabel('Y-axis');





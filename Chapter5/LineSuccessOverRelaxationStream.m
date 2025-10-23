% Leander Tenbarge: Hoffman CFD Chapter 5, Numerical Methods for Elliptic partial differential equations:
% Line Successive Over Relaxation (LSOR) for Stream Function:
% Page 175 - 176:

clear all

% System Setup as listed in page 155:
width = 5; 
height = 5; 
numPW = 100; % Number points in the width:
numPH = 100; % Number of points in the height:
dx = width/numPW;
dy = height/numPH;
n = 5000;
beta = dx/dy;
Tupper = 0;   % rankine
Tleft  = 0;   % rankine
Tright = 100; % rankine
Tlower = 100; % rankine
alpha = -2*(1+ beta^2); 
omega = 1;

% Implementing boundary conditions:
Psi = zeros(numPH,numPW); 
Ldist = linspace(0,width,numPW);
Rdist = linspace(0,height,numPH);

% Lower boundary:
for i = 1:size(Ldist,2)
    if Ldist(i) <= 1
        Psi(end,i) = 0;

    elseif Ldist(i) > 2
        Psi(end,i) = 100; 
    else
        Psi(end,i) = Ldist(i) * 100 - 100;
    end
end

% Right Boundary:
for i = 1:size(Rdist,2)
    if Rdist(i) <= 1
        Psi(i,end) = 0;

    elseif Rdist(i) > 2
        Psi(i,end) = 100; 
    else
        Psi(i,end) = Rdist(i) * 100 - 100;
    end
end


% Implementing the Line Successive Over Relaxation method to Determine Stream Function:

% Creating the A matricies:
A = spdiags([omega,alpha,omega],[-1,0,1],numPH-2,numPH-2);

% Looping over time:
for iter = 1:n
    for j = 2:numPW-1

        % Generating the Matricies:
        B = zeros(numPH-2,1);

        % Enforcing the Boundary conditions:
        B(1) = -omega*Psi(1,j);
        B(numPH-2) = -Psi(end,j)*omega;
        
        % Creating the rest of the left hand side:
        for i = 1:numPH-2
            B(i) = B(i) + alpha*(1-omega)*Psi(i,j) + (-omega*beta^2)*(Psi(i+1,j+1) + Psi(i+1,j-1));
        end

        % Solving the Tridiagonal system and replacing the column vector:
        Psi(2:numPH-1,j) = A\B;

    end
end

% Calculate velocities from stream function using the gradient function:
[dPsi_dy, dPsi_dx] = gradient(Psi, dy, dx);
u = dPsi_dy;      
v = -dPsi_dx;     
V_mag = sqrt(u.^2 + v.^2);  

% Plot velocity magnitude
figure;
imagesc(V_mag)
colorbar;
title('Velocity Magnitude');
xlabel('X Position');
ylabel('Y Position');
axis equal tight;


% Leander Tenbarge, Inviscid Burger's equation
% Anderson's CFD, page: 180
% The Rusanov (Burstein-Mirin) method:
clear all

% Initial parameters:
dx = .0025;
dt = 0.001;
xmax = 4;
nx = xmax/dx;
nt = 1000;

% Calculating the Courant number
CFL = (dt/dx);

% Setting up initial conditions:
x = linspace(0, xmax, nx);
u = zeros(nx, nt);
u(1:round(nx/2), 1) = 1;
leftValue = 1;
rightValue = 0;
omega = ((4 * CFL^2 - CFL^4) + 3)/2;

% Making the anonymous function for f(u) = u^2/2
F = @(u) u.^2/2;

% Utilizing the Rusanov (Burstein-Mirin) method
for j = 1:nt-1
    u1 = zeros(nx, 1);
    u2 = zeros(nx, 1);
    
    % Stage 1: Compute u1 at cell interfaces (j+1/2)
    for i = 1:nx-1
        u1(i) = 0.5 * (u(i+1,j) + u(i,j)) - (1/3) * (dt/dx) * (F(u(i+1,j)) - F(u(i,j)));
    end
    
    % Stage 2: Compute u2 at cell centers using u1 fluxes
    for i = 2:nx-1
        u2(i) = u(i,j) - (2/3) * (dt/dx) * (F(u1(i)) - F(u1(i-1)));
    end
    
    % Stage 3: Final update with high-order flux and dissipation
    for i = 3:nx-2
        primary = (1/24) * (dt/dx) * (-2 * F(u(i+2,j)) + 7 * F(u(i+1,j)) - 7 * F(u(i-1,j)) + 2 * F(u(i-2,j)));
        secondary = (3/8) * (dt/dx) * (F(u2(i+1)) - F(u2(i-1)));
        tertiary = (omega/24) * (u(i+2,j) - 4 * u(i+1,j) + 6 * u(i,j) - 4 * u(i-1,j) + u(i-2,j));
        u(i,j+1) = u(i,j) - primary - secondary - tertiary;
    end
    
    % Boundary conditions (simple extrapolation or hold values)
    u(1,j+1) = leftValue;
    u(2,j+1) = u(3,j+1);
    u(nx-1,j+1) = u(nx-2,j+1);
    u(nx,j+1) = rightValue;
end

% Displaying the Results
figure;
hold on
title("Rusanov (Burstein-Mirin) method, Inviscid Burgers Equation");
plot(x, u(:, 1), '-', 'DisplayName', ['t = ', num2str(dt*1, '%.3f'), ' s'], 'LineWidth', 1.5);
plot(x, u(:, round(nt/4)), '-', 'DisplayName', ['t = ', num2str(dt*round(nt/4), '%.3f'), ' s'], 'LineWidth', 1.5);
plot(x, u(:, round(nt/2)), '-', 'DisplayName', ['t = ', num2str(dt*round(nt/2), '%.3f'), ' s'], 'LineWidth', 1.5);
plot(x, u(:, round(2*nt/3)), '-', 'DisplayName', ['t = ', num2str(dt*round(2*nt/3), '%.3f'), ' s'], 'LineWidth', 1.5);
plot(x, u(:, nt), '-', 'DisplayName', ['t = ', num2str(dt*nt, '%.3f'), ' s'], 'LineWidth', 1.5);
xlabel('x');
ylabel('u');
legend('show', 'Location', 'best');
grid on
hold off
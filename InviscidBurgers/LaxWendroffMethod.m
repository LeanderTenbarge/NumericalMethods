% Leander Tenbarge, Inviscid burger's equation
% Andersons CFD, page: 179
% Lax-Wendroff method

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
x = linspace(0,xmax,nx);
u = zeros(nx,nt);
u(1:round(nx/2),1) = 1;
leftValue = 1;
rightValue = 0;

% Making the anonymous funtion for f(u) = u^2/2 and A(u) = dF(u)/du
F = @(u) u^2/2;
A = @(u1,u2) (u1+u2)/2;

% Utilizing the Lax-Wendroff Method:
for j = 1:nt-1
    for i = 2:nx-1
        u(i,j+1) = u(i,j) - (dt/(2*dx)) * (F(u(i+1,j)) - F(u(i-1,j))) + (1/2) * (dt/dx) ^ 2 * (A(u(i,j), u(i+1,j)) * F(u(i+1,j))-F(u(i,j)) - A(u(i,j),u(i-1,j))*(F(u(i,j))-F(u(i-1,j))));
    end
end

% Displaying the results:
% Displaying the Results
hold on
title("Lax-Wendroff Method, Inviscid Burgers Equation:")
plot(x,u(:, round(1)),    '-', 'DisplayName', ['t = ', num2str(dt*round(1), '%.3f'), ' s']);
plot(x, u(:, round(nt/4)),    '-', 'DisplayName', ['t = ', num2str(dt*round(nt/4), '%.3f'), ' s']);
plot(x, u(:, round(nt/2)),    '-', 'DisplayName', ['t = ', num2str(dt*round(nt/2), '%.3f'), ' s']);
plot(x, u(:, round(2*nt/3)),  '-', 'DisplayName', ['t = ', num2str(dt*round(2*nt/3), '%.3f'), ' s']);
plot(x, u(:, round(nt)),      '-', 'DisplayName', ['t = ', num2str(dt*round(nt), '%.3f'), ' s']);

xlabel('x');
ylabel('u');
legend('show', 'Location', 'best');
hold off


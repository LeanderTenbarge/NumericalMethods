% Leander Tenbarge:
% Second order, one - dimensional, wave equations:
% Page 203-204

clear all

%% Setting up inital and boundary conditions:

% Initial parameters:
dx = .5;
dt = 0.002;
xmax = 300;
nx = xmax/dx;
nt = 1000;
a = 250;
c = a*(dt/dx);

% Setting up initial conditions:
x = linspace(0,xmax,nx);
u = zeros(nx,nt);
leftValue = 0;
rightValue = 0;

% Running over the piecewise function:
for i = 1:nx
    if x(i) >= 100 && x(i) <= 220
        u(i) = 100*sin(pi*((x(i)-100)/(120)));
    end
end

% Creating the second set of inital conditions:
for i = 2:nx-1
    u(i,2) = u(i,1) + .5 * c ^ 2 * (u(i-1,1) - 2 * u(i,1) + u(i+1,1));
end


%% Implementing the Methods:
for j = 2:nt - 1
    for i = 2:nx - 1
        u(i,j+1) = 2 * u(i,j) - u(i,j-1) + c ^ 2 * (u(i-1,j) - 2 * u(i,j) + u(i+1,j));
    end
end

imagesc(u)

%% Visualization:
figure;
imagesc(x, (0:nt-1)*dt, u');
xlabel('Position x');
ylabel('Time t');
title('Wave Equation Solution - Leapfrog Method');
colorbar;

% Optional: Animate the solution
figure;
for j = 1:10:nt
    plot(x, u(:,j), "o");
    xlabel('Position x');
    ylabel('u(x,t)');
    title(sprintf('Time = %.4f', (j-1)*dt));
    ylim([-150 150]);
    grid on;
    drawnow;
    pause(0.05);
end
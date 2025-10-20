% Leander Tenbarge: Flux Corrected Transport (Predictor -> corrector) for a
% linear PDE
% Hoffman CFD: page 233 -> 234:

clear all
%% Setting up the Initial and Boundary Conditions: 

% Initial parameters:
dx = .1;
dt = 0.0002;
xmax = 400;
nx = xmax/dx;
nt = 15000;
a = 250;
c = a*(dt/dx);

% Setting up initial conditions:
x = linspace(0,xmax,nx);
u = zeros(nx,nt);
leftValue = 0;
rightValue = 0;
Uint = zeros(size(u,1),1);

% Running over the piecewise function:
for i = 1:nx
    if x(i) >= 50 && x(i) <= 110
        u(i) = 100*sin(pi*((x(i)-50)/(60)));
    end
end
uint = zeros(size(u,1),1);



%% Preparing the Flux Corrected Transport Script:

% Calculating the e1 and e2 parameters:
e1 = (1/6) * (1 + 2 * c^2);
e2 = (1/6) * (1 - c^2);

% Running over the loop:
for j = 1:nt-1
    for i = 2:nx-1
        uint(i) = u(i,j) - .5 * c * (u(i+1,j) - u(i-1,j)) + (e1 + .5 * c^2) * (u(i+1,j) - 2 * u(i,j) + u(i-1,j));
    end
    for i = 2:nx-1
        u(i,j+1) = uint(i) - e2 * (uint(i+1) - 2 * uint(i) + uint(i-1));
    end
end


% Plotting:
hold on
title("Flux Correct Transport Predictor Corrector Implementation, Courant:", num2str(c))

plot(x,u(:,2500),'o');
plot(x,u(:,5000),'o');
plot(x,u(:,7500),'o');
plot(x,u(:,10000),'o');
plot(x,u(:,12500),"o")
legend(['t = ', num2str(2500/nt) ],['t = ', num2str(5000/nt) ],['t = ', num2str(7500/nt) ],['t = ', num2str(10000/nt) ],['t = ', num2str(12500/nt) ])
hold off
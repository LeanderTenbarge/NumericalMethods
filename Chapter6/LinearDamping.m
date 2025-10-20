% Linear Hyperbolic equations:
% Application of second order viscous dampers:
% Hoffman CFD page 231 - 233:
clear all
%% Setting up the Initial and Boundary Conditions: 

% Initial parameters:
dx = .1;
dt = 0.0001;
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

% Running over the piecewise function:
for i = 1:nx
    if x(i) >= 50 && x(i) <= 110
        u(i) = 100*sin(pi*((x(i)-50)/(60)));
    end
end


% The First Plot:
Epsilon  = 0.2;
u1 = u;
for j = 1:nt-1
    for i = 2:nx-1
       u1(i,j+1) = u1(i,j) - (c/2)*(u1(i+1,j) - u1(i-1,j)) + (c^2/2)*(u1(i+1,j) - 2*u1(i,j) + u1(i-1,j)) + Epsilon*(u1(i+1,j) - 2*u1(i,j) + u1(i-1,j));
    end
end

hold on
title("Linear Damping, Lax-Wendroff Implementation for Epsilon:", num2str(Epsilon))

plot(x,u1(:,2500),'o');
plot(x,u1(:,5000),'o');
plot(x,u1(:,7500),'o');
plot(x,u1(:,10000),'o');
plot(x,u1(:,12500),"o")
legend(['t = ', num2str(2500/nt) ],['t = ', num2str(5000/nt) ],['t = ', num2str(7500/nt) ],['t = ', num2str(10000/nt) ],['t = ', num2str(12500/nt) ])
hold off


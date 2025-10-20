% Linear Hyperbolic equations:
% One dimensional Linear Examples:
% Hoffman CFD page 186 - 202:

%% Setting up the Initial and Boundary Conditions: 

% Initial parameters:
dx = .5;
dt = 0.002;
xmax = 400;
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
    if x(i) >= 50 && x(i) <= 110
        u(i) = 100*sin(pi*((x(i)-50)/(60)));
    end
end

%% Explicit methods:

% The First upwind differencing method:
fudu = u;  
for j = 1:nt-1
    for i = 2:nx
        fudu(i,j+1) = fudu(i,j) - c*(fudu(i,j) - fudu(i-1,j));
    end
end

% The Lax method:
lu = u;
for j = 1:nt-1
    for i = 2:nx-1
        lu(i,j+1) = .5 * (lu(i+1,j)+lu(i-1,j)) - (c / 2) * (lu(i+1,j) - lu(i-1,j));
    end
end

% Midpoint leapfrog method:
mlu = u;
mlu(:,2) = mlu(:,1);
for j = 2:nt-1
    for i = 2:nx-1
        mlu(i,j+1) = -c * (mlu(i+1,j)- mlu(i-1,j)) + mlu(i,j-1);
    end
end

% The Lax - Wendroff method:
lwu = u;
for j = 1:nt-1
    for i = 2:nx-1
        lwu(i,j+1) = lwu(i,j) - a*dt*((lwu(i+1,j) - lwu(i-1,j))/(2*dx)) + .5 * a^2 * dt^2 * ((lwu(i+1,j) - 2*lwu(i,j) + lwu(i-1,j))/(dx^2));
    end
end

%% Implicit methods:

% Euler BTCS:
btcsu = u;
btcsA = spdiags([.5*c, -1, -.5*c],[-1,0,1], nx, nx);
btcsB = zeros(nx,1);
for j = 1:nt-1
        
        btcsB(1) = btcsu(1,j) + leftValue;
        btcsB(end) = btcsu(end,j) + rightValue;

    for i = 2:nx-1
        btcsB(i) = -btcsu(i,j);
    end 
    btcsu(:,j+1) = btcsA\btcsB;
end

% Implicit first upwind differencing method:
udiu = u;
udiA = spdiags([-c, (1+c)], [-1, 0], nx, nx);
udiB = zeros(nx,1);
for j = 1:nt-1
    for i = 1:nx
        udiB(i) = udiu(i,j);
    end
    udiu(:,j+1) = udiA \ udiB;
end

% The Crank Nicholson method:
ucn = u;
cnA = spdiags([.25, -1,-.25], [-1,0,1], nx, nx);
cnB = zeros(nx,1);
for j = 1:nt-1
    for i = 2:nx-1
        cnB(i) = -ucn(i,j) + .25 * c * (ucn(i+1,j)-ucn(i-1,j));
    end
    cnB(1) = leftValue;
    cnB(end) = rightValue;
    ucn(:,j+1) = cnA\cnB;
end

%% Multi-step methods:

% Richtmyer/Lax-Wendroff Multi-Step Method:
rlwu = u;
uint = zeros(nx,1);
for j = 1:nt-1
    for i = 1:nx-1
        uint(i) = .5 * (rlwu(i+1,j) + rlwu(i,j)) - (.5 * c) * (rlwu(i+1,j) - rlwu(i,j));
    end
    for i = 2:nx
        rlwu(i,j+1) = rlwu(i,j) - c * (uint(i)- uint(i-1));
    end
end

% The MacCormack Method:
mcu = u;
ustar = zeros(nx,1);
for j = 1:nt-1
    % Predictor step:
    for i = 2:nx-1
        ustar(i) = mcu(i,j) - c * (mcu(i+1,j)-mcu(i,j));
    end
    % Corrector step:
    for i = 2:nx-1
        mcu(i,j+1) = .5 * ((mcu(i,j)+ustar(i)) - c * (ustar(i) - ustar(i-1)));
    end
end


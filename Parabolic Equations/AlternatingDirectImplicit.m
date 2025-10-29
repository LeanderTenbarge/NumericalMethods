% Leander Tenbarge
% Alternating Direct Implicit, Approximate Factorization Method for 2d unsteady heat equation:
% Page 137 - 141: Anderson's CFD


%% Initial Parameters:
alpha = 1.7917e-04;
dx = 0.01;    
dy = 0.01;    
dt = 1;    
Nx = 100;     
Ny = 100;     
Nt = 1000;     
To = 0.00;    
T1 = 0.00;    % Temperature Top 
T2 = 0.00;    % Temperature Right
T3 = 300.0;    % Temperature Bottom 
T4 = 0.00;    % Temperature Left
T = To*ones(Nx,Ny,Nt); 
Residuals = zeros(Nt,1);

% Enforcing Inital and boundary conditions:
T(1,:,:) = T1;
T(:,end,:) = T2;
T(:,1,:) = T4;
T(end,:,:) = T3;
Tstar = zeros(Nx-2,Ny-2); 



%% Solution:

% Initializing the solution parameters:
rx = alpha * (dt)/(dx^2);
ry = alpha * (dt)/(dy^2);

% Creating the Sparse Diagonals (Only For Interior Nodes):
A1 = spdiags([(-rx/2),(1+rx),(-rx/2)],[-1,0,1],Nx-2,Nx-2); % Length of columns is Nx:
A2 = spdiags([(-ry/2),(1+ry),(-ry/2)],[-1,0,1],Ny-2,Ny-2); % Length of columns is Ny:

% The Solution loop:
for t = 1:Nt-1

    % Conducting the X sweep (Explicit over the jth column):
    for j = 2:Ny-1

        % Determining the B array(offset in i by i-1):
        b = zeros(Nx-2,1);
        for i = 1:Nx-2
            b(i) = b(i) + rx * (T(i,j,t) -2 * T(i+1,j,t) + T(i+2,j,t)) + ry * (T(i+1,j-1,t) - 2 * T(i+1,j,t) + T(i+1,j+1,t));
        end

        % Determining the Solution:
        Tstar(:,j-1) = A1\b;
    end

    % Conducting the Y sweep (Explicit over the ith column)
    for i = 2:Nx-1

        % Determining the B array(offset in j by j-1):
        b = zeros(Ny-2,1);
        for j = 1:Ny-2
            b(j) = b(j) + Tstar(i-1,j);
        end

        % Determining the Solution:
        T(i,2:Ny-1,t+1) =  A2\b + T(i,2:Ny-1,t)';
    end
    Residuals(t) = max(abs(T(:,:,t+1) - T(:,:,t)), [], 'all');
end



%% Compute Analytical Solution

% Centerline position
x_center = L/2;
y_vals = linspace(0, W, Ny);

% Initialize centerline temperature
T_centerline = zeros(Ny, 1);

% Compute analytical solution at centerline
for j = 1:Ny
    y = y_vals(j);
    
    % Sum the series
    for n = 1:Ny
        sinh_num = sinh(n*pi*y/L);
        sinh_den = sinh(n*pi*W/L);
        T_centerline(j) = T_centerline(j) + T3 * 2.0 * (1 - cos(n*pi))/(n*pi) * (sinh_num/sinh_den) * sin(n*pi*x_center/L);
    end
end


Tmin = min(T(:));
Tmax = max(T(:));

%Subplot(1: t = 50*dt)
subplot(2,3,1);
imagesc(T(:,:,50), [Tmin Tmax]);
colorbar;
title("T = " + 50*dt + " sec");
xlabel('X Index');
ylabel('Y Index');
axis equal tight;

% Subplot 2: t = 100*dt
subplot(2,3,2);
imagesc(T(:,:,100), [Tmin Tmax]);
colorbar;
title("T = " + 100*dt + " sec");
xlabel('X Index');
ylabel('Y Index');
axis equal tight;


% Subplot 3: t = Nt*dt (final time)
subplot(2,3,3);
imagesc(T(:,:,end), [Tmin Tmax]);
colorbar;
title("T = " + Nt*dt + " sec");
xlabel('X Index');
ylabel('Y Index');
axis equal tight;

% Subplot 4: Residuals:
subplot(2,3,4)
semilogy(1:Nt,Residuals)
title("Residuals")
xlabel('Iterations');
ylabel('Maximum Residual');
axis tight;

% Subplot 5 (Centerline Temperature dist
subplot(2,3,5)
hold on
plot(y_vals,T_centerline);
plot(y_vals,T(:,round(Ny/2),end));
title('Steady State Temperature Distribution (k)')
xlabel('X dimension (ith)');
ylabel('Temperature');
% Overall figure title
sgtitle("Temperature Distribution at Various Time Steps");



clc; clear; close;

clear all
% Parameters
alpha = 0.000217;        % Kinematic viscosity (m^2/s)
U0 = 40;                 % Velocity of lower plate (m/s)
h = 0.04;                % Distance between plates (m)
Delta_x = 0.001;         % Spatial step (m)
x = 0:Delta_x:h;         % y grid points
n_x = length(x);         % Number of grid points in y

% Time parameters
Delta_t1 = 0.002;        % Time step 1 (s)
Delta_t2 = 0.00232;      % Time step 2 (s)
t_end = 10;              % Final time (s)

% Number of time steps for each case
n_t1 = round(t_end / Delta_t1);
n_t2 = round(t_end / Delta_t2);

% Initialize velocity
u_init = zeros(n_x, 1);
u_init(1) = U0;          % Boundary condition at y = 0

% Function for Crank-Nicolson scheme
function u = crank_nicolson(u, n_t, dt, dx, alpha, U0)
    r = alpha * dt / dx^2;
    N = length(u);

    % Construct A (implicit) and B (explicit) matrices for interior nodes
    main_diag_A = (1 + r) * ones(N-2, 1);
    off_diag_A = -0.5 * r * ones(N-3, 1);
    A = diag(main_diag_A) + diag(off_diag_A, 1) + diag(off_diag_A, -1);

    main_diag_B = (1 - r) * ones(N-2, 1);
    off_diag_B = 0.5 * r * ones(N-3, 1);
    B = diag(main_diag_B) + diag(off_diag_B, 1) + diag(off_diag_B, -1);

    for t = 1:n_t
        % RHS = B * u_old_interior + boundary terms
        u_interior = u(2:end-1);
        rhs = B * u_interior;

        % Add boundary terms
        rhs(1) = rhs(1) + 0.5 * r * U0; % Left boundary
        rhs(end) = rhs(end) + 0;        % Right boundary (u=0)

        % Solve linear system
        u_new_interior = A \ rhs;

        % Update solution including boundaries
        u = [U0; u_new_interior; 0];
    end
end

% Run simulations for both Δt
u1 = crank_nicolson(u_init, n_t1, Delta_t1, Delta_x, alpha, U0);
u2 = crank_nicolson(u_init, n_t2, Delta_t2, Delta_x, alpha, U0);

% Plot results
figure;
plot(x, u1, '-o', 'DisplayName', 'Crank-Nicolson, Δt = 0.002 s');
hold on;
plot(x, u2, '-x', 'DisplayName', 'Crank-Nicolson, Δt = 0.00232 s');
xlabel('y (m)');
ylabel('Velocity u (m/s)');
legend show;
title('Velocity Profile between Parallel Plates (Crank-Nicolson)');
grid on;

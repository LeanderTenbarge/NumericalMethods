% Leander Tenbarge, Approximate Factorized ADI method 
% Method as layed out in Hoffman CFD page no.86:

clear all

%% Initial Parameters:

% Constants and such:
dx = 0.001;    % ft
dy = 0.001;    % ft
dt = 0.01;    % sec
Nx = 100;     % Number of X nodes
Ny = 100;     % Number of Y nodes
Nt = 500;     % Number of Temperature nodes
To = 0;    % Inital Condition
T1 = 200.0;  % Temperature left 
T2 = 200.0;  % Temperature lower
T3 = 0;    % Temperature Right 
T4 = 0;    % Temperature Top
alpha = 1.79e-04; 


% Creating the Domain & enforcing boundary conditions:
T = zeros(Nx,Ny,Nt);
T(1,:,:) = T1;
T(:,end,:) = T2;
T(end,:,:) = T3;
T(:,1,:) = T4;

%% Preliminary calculations
% Calculating the Variables (Refer to sourcematerial):
deltax = (alpha*dt)/(2*dx^2); 
deltay = (alpha*dt)/(2*dy^2);

% Creating the x-sweep solution matricies:
A1 = spdiags([(-deltax/2), (1+deltax), (-deltax/2) ],[-1, 0, 1 ], size(T,2)-2, size(T,2)-2);

% Creating the y-sweep solution matricies:
A2 = spdiags([(-deltay/2),(1+deltay),(-deltay/2)], [-1, 0, 1 ], size(T,1)-2, size(T,1)-2);

%The temporary solution arrays:
B1 = zeros(size(T,2)-2,1);
B2 = zeros(size(T,1)-2,1);



%% The Solution loop over time:
for n = 1:size(T,3)-1    
    % Creating the x-sweep:

    % Creating the intermediate temperature field:
    Tint = T(:,:,n);
    
    % Initializing the temporary B1 array
    B1 = zeros(size(T,2)-2,1);
    
    
    % Performing the x-sweep
    for j = 2:size(T,2)-1
    
        % Enforcing edge conditions before the loop:
    
        B1(1) = (deltax/2)*Tint(1,j);
        B1(end) = (deltax/2)*Tint(end,j);
    
        for i = 2:size(T,1)-1
    
            % Calculating the B1 array for each tridiagional system 
            B1(i-1) = B1(i-1)+ (1-deltay)*T(i,j,n) + (deltay/2)*(T(i,j+1)+T(i,j-1));
    
        end
        
        % Calculating the results and inserting them into the temporary array:
        Tint(2:size(T,2)-1,j) = (A1\B1);
    
    end
    
    
    % Creating the y-sweep:
    
    % Initializing the inital B2 array:
    B2 = zeros(size(T,1)-2,1);
    
    % Performing the y-sweep:
    for i = 2:size(T,1)-1
    
        % Enforcing edge conditions before the loop:
        B2(1) = (deltay/2)*Tint(i,1);
        B2(end) = (deltay/2)*Tint(i,end);
    
        for j = 2:size(T,2)-1
    
            % Calculating the B2 array for each tridiagional system 
            B2(j-1) = (1-deltax)*Tint(i,j)+ (deltax/2)*(Tint(i+1,j)+Tint(i-1,j));
    
        end
        % Calculating the results and inserting them into the system:
        T(i,2:size(T,2)-1,n+1) = (A2\B2);
    
    end
end

imagesc(T(:,:,50))
% Leander Tenbarge, Approximate Factorized ADI method 
% Method as layed out in Hoffman CFD page no.86:

clear all

% Initial Parameters:
dx = 0.1;    % ft
dy = 0.1;    % ft
dt = 0.1;    % sec
Nx = 35;     % Number of X nodes
Ny = 35;     % Number of Y nodes
Nt = 50;     % Number of Temperature nodes
To = 0.0;    % Inital Condition
T1 = 200.0;  % Temperature left 
T2 = 200.0;  % Temperature lower
T3 = 0.0;    % Temperature Right 
T4 = 0.0;    % Temperature Top
time = 50;    % sec
T = zeros(Nx,Ny,Nt); % The temperature matrix
alpha = 0.645/3600; % ft^2/hr -> ft^2/sec

% Enforcing Inital and boundary conditions:
T(1,:,:) = T1;
T(:,end,:) = T2;
T(end,:,:) = T3;
T(:,1,:) = T4;

%Calculating the parameters
alpha = 0.645/3600; % ft^2/hr -> ft^2/sec
Dx = (alpha*dt)/(dx^2);
Dy = (alpha*dt)/(dy^2);

% Constructing the A & B matricies and fixing the column offset issue:
b1 = zeros(size(T,1)-2, size(T,2)-2);
A1 = spdiags([-Dx/2,(1+Dx),-Dx/2], [-1,0,1], numel(b1),numel(b1));
b2 = zeros(size(T,1)-2, size(T,2)-2);
A2 = spdiags([-Dy/2,(1+Dy),-Dy/2], [-1,0,1], numel(b1),numel(b1));

% Routine for fixing column looping issues:
for column = 1:size(b1,2)
    base = (column - 1)*size(b1,2);
    if column < size(b1,2)
        A1(base+size(b1,1) + 1, base+size(b1,1)) = 0;
        A1(base+size(b1,1), base+size(b1,1)+ 1 ) = 0;
        A2(base+size(b1,1) + 1, base+size(b1,1)) = 0;
        A2(base+size(b1,1), base+size(b1,1)+ 1 ) = 0;
    end
end



% Looping over the matrix creation and solving;
for n = 1:size(T,3)-1

% Constructing the b1 matricies
    for i = 1:size(b1,1)
        for j = 1:size(b1,2)
            % Top boundary (row above current point)
            if i == 1
                b1(i,j) = b1(i,j) + (Dx/2)*(T(i,j+1,n));
            end
            % Bottom boundary (row below current point)
            if i == size(b1,1)
                b1(i,j) = b1(i,j) + (Dx/2)*(T(i+2,j+1,n));
            end
            
            % Constructing the rest of the elements
            b1(i,j) = b1(i,j) + (1-Dy)*T(i+1,j+1,n) + (Dy/2)*(T(i+1,j+2,n))+(Dy/2)*(T(i+1,j));
        end
    end
    

    % Flattening the B Matricies:
    b1 = reshape(b1, [], 1);
    x1 = A1\b1;
    x1 = reshape(x1,size(T,1)-2, size(T,2)-2);
    
    % Creating the temperature field at n+=1/2 by embedding the x1:
    TempT = T(:,:,n);
    TempT(2:34,2:34) = x1;

    % Contructing the b2 matricies:
    for i = 1:size(b2,1)
        for j = 1:size(b2,2)
            % Left boundary (column next to current point)
            if j == 1
                b2(i,j) = b2(i,j) + (Dy/2)*(T(i+1,j,n));
            end

            % Right Boundary (column next to current point)
            if j == size(b1,2)
                b2(i,j) = b2(i,j) + (Dy/2)*(T(i+1,j+2,n));
            end
            
            % Constructing the rest of the elements (in terms of)
            b2(i,j) = b2(i,j) + (1-Dx)*TempT(i+1,j+1) + (Dx/2)*(TempT(i+2,j+1))+(Dx/2)*(TempT(i,j+1));
            
        end
    end
    
    b2 = reshape(b2, [], 1);
    x2 = A2\b2;
    x2 = reshape(x2,size(T,1)-2, size(T,2)-2);
    
    T(2:size(x2,1)+1, 2:size(x2,1)+1 ,n+1) = x2;


    % Reset & Unflatten the matricies
    b1 = reshape(b1,size(T,1)-2, size(T,2)-2);
    b2 = reshape(b2,size(T,1)-2, size(T,2)-2);
end


imagesc(T(:,:,n))
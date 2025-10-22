% Leander Tenbarge, Hoffman CFD Chapter 5, Numerical Methods for Elliptic partial differential equations:
% Implicit 5 point formulation:
% Page 153-155

clear all

% System Setup as listed in page 155:
width = 1;  % feet
height = 1; % feet
numPointsWidth = 150;
numPointsHeight = 150;
dx = width/numPointsWidth;
dy = height/numPointsHeight;
beta = dx/dy;
Tupper = 0; % rankine
Tleft  = 0; % rankine
Tright = 100; % rankine
Tlower = 100; % rankine


% Basic 5 grid method (Pentadiagonal Matrices, Non-iterative formulation):
T = zeros(numPointsHeight,numPointsWidth); % Generating the solution space
T(:,1) = Tleft;
T(1,:) = Tupper;
T(:,end) = Tright;
T(end,:) = Tlower;
alpha = -2*(1+ beta^2);


% Constructing the b matrix in Ax = b:
b = zeros(size(T,1)-2, size(T,2)-2);  % (ny-2)x(nx-2)

for i = 1:size(b,1)
    for j = 1:size(b,2)

        % Top boundary (row above current point)
        if i == 1
            b(i,j) = b(i,j) - T(i,j+1);
        end

        % Bottom boundary (row below current point)
        if i == size(b,1)
            b(i,j) = b(i,j) - T(i+2,j);
        end

        % Left boundary (column to the left of current point)
         if j == 1 
            b(i,j) = b(i,j) - (beta^2)*T(i+1,j);
         end

        % Right boundary (column to the right of current point)
         if j == size(b,2)
            b(i,j) = b(i,j) - (beta^2)*T(i+1,j+2);
        end

    end
end

% Creating the A matricies, via the Spdiags function and correcting for the Column jump effects:
% Because we need to adjust for the columns wrapping along the beta terms are a length(column) number of points laterally,
% as represented in the matrix.

A = spdiags([beta^2, 1, alpha, 1, beta^2],[-size(b,1),-1,0,1,size(b,1)],numel(b),numel(b)); % Contructing the sparse matricies:

% Cleaning up the vertical wrap around:
for column = 1:size(b,2)
    base = (column - 1)*size(b,1);
    if column < size(b,2)
        % Remove wrap-around entries between adjacent columns
        A(base + size(b,1), base + size(b,1) + 1) = 0;
        A(base + size(b,1) + 1, base + size(b,1)) = 0;
    end
end

% Flattening the B Matricies
b = reshape(b, [], 1);

% Solving the linear equation,and reformatting back into correct dimensions
x = A \ b;
x = reshape(x,numPointsHeight-2,numPointsWidth-2);

% Plot heat map
imagesc(x);          % display the matrix as an image
colorbar;               % show color scale
colormap('jet');        % choose a colormap ('jet', 'parula', 'hot', etc.)
title('Five-Point direct formulation for Steady State Heat equation');
xlabel('X-axis');
ylabel('Y-axis');


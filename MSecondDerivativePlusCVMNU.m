function laplacian = MSecondDerivativePlusCVMNU(n, delta, boundaryConditionLeft, boundaryConditionRight)
    % Computes the second derivative operator matrix on a non-uniform grid
    % with Dirichlet (0) or Neumann (1) boundary conditions.
    % delta: vector of spacing between points (length n-1)

    % Step 1: Reconstruct the grid
    x = cumsum([0; delta(:)]);  % Grid points
    laplacian = zeros(n);

    % Step 2: Fill interior rows using finite differences for non-uniform grids
    for i = 2:n-1
        h_i = x(i) - x(i-1);
        h_ip1 = x(i+1) - x(i);
        denom = h_i * h_ip1 * (h_i + h_ip1);

        laplacian(i, i-1) = 2 / (h_i * (h_i + h_ip1));
        laplacian(i, i)   = -2 / (h_i * h_ip1);
        laplacian(i, i+1) = 2 / (h_ip1 * (h_i + h_ip1));
    end

    % Step 3: Handle boundaries
    if boundaryConditionLeft == 0  % Dirichlet
        laplacian(1,:) = 0;
    elseif boundaryConditionLeft == 1  % Neumann (forward second difference)
        h = x(2) - x(1);
        laplacian(1,1) = 1 / h^2;
        laplacian(1,2) = -2 / h^2;
        laplacian(1,3) = 1 / h^2;
    end

    if boundaryConditionRight == 0  % Dirichlet
        laplacian(end,:) = 0;
    elseif boundaryConditionRight == 1  % Neumann (backward second difference)
        h = x(end) - x(end-1);
        laplacian(end,end)   = 1 / h^2;
        laplacian(end,end-1) = -2 / h^2;
        laplacian(end,end-2) = 1 / h^2;
    end
end







% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% function laplacian = MSecondDerivativePlusCVMNU(n, delta, boundaryConditionLeft, boundaryConditionRight)
%     % MSecondDerivativePlusC computes the second derivative of the matrix L
%     % with specified boundary conditions and a boundary matrix.
%     % Arguments:
%     % L - the input matrix
%     % delta - the spacing (step size h)
%     % boundaryConditionLeft - 0 (Dirichlet) or 1 (Neumann) for the left boundary
%     % boundaryConditionRight - 0 (Dirichlet) or 1 (Neumann) for the right boundary
%     % boundaryMatrix - matrix of boundary correction values
% 
%     % [n, ~] = size(L);  % size of the matrix
%     h = delta;         % step size
% 
%     % Create the Laplacian matrix for second derivative (central difference)
%     laplacian = -2 * diag(ones(n, 1)) + diag(ones(n-1, 1), 1) + diag(ones(n-1, 1), -1);
% 
%     % Apply boundary conditions for left boundary
%     if boundaryConditionLeft == 0  % Dirichlet boundary condition
%         %laplacian(1, 1) = 1;
%         %laplacian(1, 2) = 0;  % Dirichlet: boundary is set to a specific value
%     elseif boundaryConditionLeft == 1  % Neumann boundary condition
%         laplacian(1, 1) = -1;
%         %laplacian(1, 2) = 2;  % Neumann: derivative at the boundary
%     end
% 
%     % Apply boundary conditions for right boundary
%     if boundaryConditionRight == 0  % Dirichlet boundary condition
%         %laplacian(n, n) = 1;
%         %laplacian(n, n-1) = 0;  % Dirichlet: boundary is set to a specific value
%     elseif boundaryConditionRight == 1  % Neumann boundary condition
%         laplacian(n, n) = -1;
%         %laplacian(n, n-1) = 2;  % Neumann: derivative at the boundary
%     end
% 
%     % Scale the Laplacian by 1/h^2 for the second derivative
%     laplacian = laplacian / (h^2);
% 
%     % % Compute the second derivative
%     % D = laplacian * L; %+ boundaryMatrix; 
%     % 
%     % % Apply the boundaryMatrix correction to boundary rows (first and last row)
%     % D(1, :) = D(1, :) + (boundaryMatrix(1, :));%/ (h^2);  % Top row correction
%     % D(end, :) = D(end, :) + (boundaryMatrix(end, :));%/ (h^2);  % Bottom row correction
% 
% end

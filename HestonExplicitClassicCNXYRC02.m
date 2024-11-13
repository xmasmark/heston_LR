function U = HestonExplicitClassicCNXYRC02(params,K,r,q,S,V,T,mode)

    %mode is to decide the system resolution:
    %0--> Euler
    %1--> Crank Nicolson with the \
    %2--> Crank Nicolson with GMRS

    % Finite differences for the Heston PDE for a European Call
    % Uses even grid sizes
    % In 'T Hout and Foulon "ADI Finite Difference Schemes for Option Pricing
    % in the Heston Modelo with Correlation" Int J of Num Analysis and Modeling, 2010.
    % Thesis by Sensi Li and paper by Vassilis Galiotos
    % INPUTS
    %    params = 6x1 vector of Heston parameters
    %    K = Strike price
    %    r = risk free rate
    %    q = Dividend yield
    %    S = vector for stock price grid
    %    V = vector for volatility grid
    %    T = vector for maturity grid
    % OUTPUT
    %    U = U(S,v) 2-D array of size (nS+1)x(nV+1) for the call price
    
    % Heston parameters
    kappa = params(1);
    theta = params(2);
    sigma = params(3);
    v0    = params(4);
    rho   = params(5);
    lambda = params(6);
    epsilon = 0.01;
    
    % Length of stock price, volatility, and maturity
    NS = length(S);
    NV = length(V);
    NT = length(T);
    Smin = S(1);  Smax = S(NS);
    Vmin = V(1);  Vmax = V(NV);
    Tmin = T(1);  Tmax = T(NT);
    
    % Increment for Stock Price, Volatility, and Maturity
    ds = (Smax-Smin)/(NS-1);
    dv = (Vmax-Vmin)/(NV-1);
    dt = (Tmax-Tmin)/(NT-1);
    
    % Initialize the 2-D grid with zeros
    U = zeros(NS,NV);
    
    % Temporary grid for previous time steps
    u = zeros(NS,NV);
    
    % Solve the PDE
    % Round each value of U(S,v,t) at each step
    % Boundary condition for t = Maturity
    
    %single loop here producing X and Y
    for s=1:NS
	    for v=1:NV
		    U(s,v) = max(S(s) - K, 0);
	    end
    end
    
    X=zeros(NS,1);
    for s=1:NS
        X(s)=max(S(s)-K,0);
    end

    Y=zeros(NV,1);
    for v=1:NV
	    Y(v)=1;
    end
    
    SMatrix = diag(S);
    S2Matrix = diag(S.^2);
    VMatrix = diag(V);

    for t = 1:NT-1

        x = [X];
        y = [Y];

        [AX,AY] = HestonMatVec(x,y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);
        [BX,BY] = HestonMatVecBoundaries(NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, K, Tmax, t, T);

        % Set GMRES parameters
        restart = 80;  % Restart after 20 iterations (example value)
        tol = 1e-5;  % Tolerance for convergence
        max_iter = 100;  % Maximum number of iterations


        % BX = [(1+r*dt/2)*x, (-dt/2)*AX,b1x,b2x,b3x,b4x,b5x];
        % BY = [           y,         AY,b1y,b2y,b3y,b4y,b5y];

        FX = [(1+r*dt/2)*x, (-dt/2)*AX,BX];
        FY = [           y,         AY,BY];

        [BXc,BYc]=CompressData(FX, FY, epsilon);

        [X, Y] = GMRES_XYv01(BXc, BYc, x, y, restart, tol, max_iter);
    end    
    U=X*Y';
end

% % % function R = SVDRank(U, epsilon)
% % % 
% % %     %builtin_rank = rank(U,epsilon);
% % % 
% % %     [uc,sc,vs]=svd(U);
% % % 
% % %     sigSum = 0;
% % % 
% % %     %the look should go to the minimum of the two dimensions
% % %     %so a step is missing and it should be added when moving
% % %     %this to a separate function
% % %     for i = 1: size(sc,2)
% % %         sigSum = sigSum + power(sc(i,i),2);
% % %     end
% % % 
% % %     threshold = power(epsilon,2)*sigSum;
% % % 
% % %     rank_home_made = 0;
% % % 
% % %     errorSquare = 0;
% % % 
% % %     for i=size(sc,2):-1:1
% % %         errorSquare = errorSquare + power(sc(i,i),2);
% % %         if errorSquare > threshold
% % %             rank_home_made = i + 1;
% % %             break;
% % %         end
% % %     end
% % % 
% % %     R=rank_home_made;
% % % end

function D = MDerivativeC(matrix, delta, boundaryConditionLeft, boundaryConditionRight, boundaryMatrix)
    % MDerivativeP computes the first derivative of the input matrix
    % with specified boundary conditions and a boundary matrix.
    % Arguments:
    % matrix - the input matrix
    % delta - the step size (h)
    % boundaryConditionLeft - 0 (Dirichlet) or 1 (Neumann) for the left boundary
    % boundaryConditionRight - 0 (Dirichlet) or 1 (Neumann) for the right boundary
    % boundaryMatrix - the boundary matrix for correction

    [records, fields] = size(matrix); % size of the matrix
    h = delta; % step size

    % Create the "Laplacian" matrix (finite difference approximation)
    temporary = -diag(ones(records-1, 1), -1) + diag(ones(records-1, 1), 1);
    
    % Handle boundary conditions
    if boundaryConditionLeft == 0  % Dirichlet left
        temporary(1,1) = 0;
    elseif boundaryConditionLeft == 1  % Neumann left
        temporary(1,1) = -1;
    end

    if boundaryConditionRight == 0  % Dirichlet right
        temporary(records, records) = 0;
    elseif boundaryConditionRight == 1  % Neumann right
        temporary(records, records) = 1;
    end

    % Scale by 1/(2*h) for central difference
    temporary = temporary / (2 * h);

    % Compute the derivative column by column
    D = temporary * matrix; % + boundaryMatrix; 

    % Apply the boundaryMatrix correction to boundary rows (first and last row)
    D(1, :) = D(1, :) + (boundaryMatrix(1, :));%/(2*h); % Top row
    D(end, :) = D(end, :) + (boundaryMatrix(end, :));%/(2*h); % Bottom row
end

function D = MDerivativeV(matrix, delta, boundaryConditionLeft, boundaryConditionRight)
    % MDerivativeP computes the first derivative of the input matrix
    % with specified boundary conditions and a boundary matrix.
    % Arguments:
    % matrix - the input matrix
    % delta - the step size (h)
    % boundaryConditionLeft - 0 (Dirichlet) or 1 (Neumann) for the left boundary
    % boundaryConditionRight - 0 (Dirichlet) or 1 (Neumann) for the right boundary
    % boundaryMatrix - the boundary matrix for correction

    [records, fields] = size(matrix); % size of the matrix
    h = delta; % step size

    % Create the "Laplacian" matrix (finite difference approximation)
    temporary = -diag(ones(records-1, 1), -1) + diag(ones(records-1, 1), 1);
    
    % Handle boundary conditions
    if boundaryConditionLeft == 0  % Dirichlet left
        temporary(1,1) = 0;
    elseif boundaryConditionLeft == 1  % Neumann left
        temporary(1,1) = -1;
    end

    if boundaryConditionRight == 0  % Dirichlet right
        temporary(records, records) = 0;
    elseif boundaryConditionRight == 1  % Neumann right
        temporary(records, records) = 1;
    end

    % Scale by 1/(2*h) for central difference
    temporary = temporary / (2 * h);

    % Compute the derivative column by column
    D = temporary * matrix; % + boundaryMatrix; 

    % Apply the boundaryMatrix correction to boundary rows (first and last row)
    % D(1, :) = D(1, :) + (boundaryMatrix(1, :));%/(2*h); % Top row
    % D(end, :) = D(end, :) + (boundaryMatrix(end, :));%/(2*h); % Bottom row
end

function temporary = MDerivativeVM(matrix, delta, boundaryConditionLeft, boundaryConditionRight)
    % MDerivativeP computes the first derivative of the input matrix
    % with specified boundary conditions and a boundary matrix.
    % Arguments:
    % matrix - the input matrix
    % delta - the step size (h)
    % boundaryConditionLeft - 0 (Dirichlet) or 1 (Neumann) for the left boundary
    % boundaryConditionRight - 0 (Dirichlet) or 1 (Neumann) for the right boundary
    % boundaryMatrix - the boundary matrix for correction

    [records, fields] = size(matrix); % size of the matrix
    h = delta; % step size

    % Create the "Laplacian" matrix (finite difference approximation)
    temporary = -diag(ones(records-1, 1), -1) + diag(ones(records-1, 1), 1);
    
    % Handle boundary conditions
    if boundaryConditionLeft == 0  % Dirichlet left
        temporary(1,1) = 0;
    elseif boundaryConditionLeft == 1  % Neumann left
        temporary(1,1) = -1;
    end

    if boundaryConditionRight == 0  % Dirichlet right
        temporary(records, records) = 0;
    elseif boundaryConditionRight == 1  % Neumann right
        temporary(records, records) = 1;
    end

    % Scale by 1/(2*h) for central difference
    temporary = temporary / (2 * h);

    % Compute the derivative column by column
    % D = temporary * matrix; % + boundaryMatrix; 

    % Apply the boundaryMatrix correction to boundary rows (first and last row)
    % D(1, :) = D(1, :) + (boundaryMatrix(1, :));%/(2*h); % Top row
    % D(end, :) = D(end, :) + (boundaryMatrix(end, :));%/(2*h); % Bottom row
end



function D = MSecondDerivativePlusC(L, delta, boundaryConditionLeft, boundaryConditionRight, boundaryMatrix)
    % MSecondDerivativePlusC computes the second derivative of the matrix L
    % with specified boundary conditions and a boundary matrix.
    % Arguments:
    % L - the input matrix
    % delta - the spacing (step size h)
    % boundaryConditionLeft - 0 (Dirichlet) or 1 (Neumann) for the left boundary
    % boundaryConditionRight - 0 (Dirichlet) or 1 (Neumann) for the right boundary
    % boundaryMatrix - matrix of boundary correction values

    [n, ~] = size(L);  % size of the matrix
    h = delta;         % step size

    % Create the Laplacian matrix for second derivative (central difference)
    laplacian = -2 * diag(ones(n, 1)) + diag(ones(n-1, 1), 1) + diag(ones(n-1, 1), -1);
    
    % Apply boundary conditions for left boundary
    if boundaryConditionLeft == 0  % Dirichlet boundary condition
        laplacian(1, 1) = 1;
        laplacian(1, 2) = 0;  % Dirichlet: boundary is set to a specific value
    elseif boundaryConditionLeft == 1  % Neumann boundary condition
        laplacian(1, 1) = -2;
        laplacian(1, 2) = 2;  % Neumann: derivative at the boundary
    end

    % Apply boundary conditions for right boundary
    if boundaryConditionRight == 0  % Dirichlet boundary condition
        laplacian(n, n) = 1;
        laplacian(n, n-1) = 0;  % Dirichlet: boundary is set to a specific value
    elseif boundaryConditionRight == 1  % Neumann boundary condition
        laplacian(n, n) = -2;
        laplacian(n, n-1) = 2;  % Neumann: derivative at the boundary
    end

    % Scale the Laplacian by 1/h^2 for the second derivative
    laplacian = laplacian / (h^2);

    % Compute the second derivative
    D = laplacian * L; %+ boundaryMatrix; 
    
    % Apply the boundaryMatrix correction to boundary rows (first and last row)
    D(1, :) = D(1, :) + (boundaryMatrix(1, :));%/ (h^2);  % Top row correction
    D(end, :) = D(end, :) + (boundaryMatrix(end, :));%/ (h^2);  % Bottom row correction

end


function laplacian = MSecondDerivativePlusCVM(L, delta, boundaryConditionLeft, boundaryConditionRight)
    % MSecondDerivativePlusC computes the second derivative of the matrix L
    % with specified boundary conditions and a boundary matrix.
    % Arguments:
    % L - the input matrix
    % delta - the spacing (step size h)
    % boundaryConditionLeft - 0 (Dirichlet) or 1 (Neumann) for the left boundary
    % boundaryConditionRight - 0 (Dirichlet) or 1 (Neumann) for the right boundary
    % boundaryMatrix - matrix of boundary correction values

    [n, ~] = size(L);  % size of the matrix
    h = delta;         % step size

    % Create the Laplacian matrix for second derivative (central difference)
    laplacian = -2 * diag(ones(n, 1)) + diag(ones(n-1, 1), 1) + diag(ones(n-1, 1), -1);
    
    % Apply boundary conditions for left boundary
    if boundaryConditionLeft == 0  % Dirichlet boundary condition
        %laplacian(1, 1) = 1;
        %laplacian(1, 2) = 0;  % Dirichlet: boundary is set to a specific value
    elseif boundaryConditionLeft == 1  % Neumann boundary condition
        laplacian(1, 1) = -1;
        %laplacian(1, 2) = 2;  % Neumann: derivative at the boundary
    end

    % Apply boundary conditions for right boundary
    if boundaryConditionRight == 0  % Dirichlet boundary condition
        %laplacian(n, n) = 1;
        %laplacian(n, n-1) = 0;  % Dirichlet: boundary is set to a specific value
    elseif boundaryConditionRight == 1  % Neumann boundary condition
        laplacian(n, n) = -1;
        %laplacian(n, n-1) = 2;  % Neumann: derivative at the boundary
    end

    % Scale the Laplacian by 1/h^2 for the second derivative
    laplacian = laplacian / (h^2);

    % % Compute the second derivative
    % D = laplacian * L; %+ boundaryMatrix; 
    % 
    % % Apply the boundaryMatrix correction to boundary rows (first and last row)
    % D(1, :) = D(1, :) + (boundaryMatrix(1, :));%/ (h^2);  % Top row correction
    % D(end, :) = D(end, :) + (boundaryMatrix(end, :));%/ (h^2);  % Bottom row correction

end


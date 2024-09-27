function U = HestonExplicitClassicCN_GMRS(params,K,r,q,S,V,T)

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
    
    %vectorised U = X*Y' where apostrophe is 
    %the sign for transposition (column to row).
    
    % test = U-X*Y';
    % 
    % matrix_rank=[];
    % matrix_rank_home_made=[];

    SMatrix = diag(S);
    S2Matrix = diag(S.^2);
    VMatrix = diag(V);

    %%%%%%%%% this works beautifully -- don't lose it!
    % % % % % for t = 1:NT-1
    % % % % %     % Initialize u for this time step
    % % % % %     u = U;
    % % % % % 
    % % % % %     % Set up boundary matrices and derivatives in S
    % % % % %     boundaryMatrix_S = zeros(NS, NV);
    % % % % %     boundaryMatrix_S(NS,:) = (1 / ds);  %Neumann condition at large S
    % % % % % 
    % % % % %     firstSDerivative = MDerivativeC(u, ds, 0, 1, boundaryMatrix_S);
    % % % % %     secondSDerivative = MSecondDerivativePlusC(u, ds, 0, 1, boundaryMatrix_S);
    % % % % % 
    % % % % %     % Set up boundary matrices and derivatives in V
    % % % % %     boundaryMatrix_V = zeros(NS, NV);
    % % % % %     boundaryMatrixT = boundaryMatrix_V';  % Transpose for volatility derivatives
    % % % % %     boundaryMatrixT(NV,:) = S;  % Boundary at max volatility
    % % % % % 
    % % % % %     % Boundary Condition when V = 0 (volatility = 0)
    % % % % %     discountedPayoff = max((S - K * exp(-r * (Tmax - T(t)))), 0);
    % % % % %     %discountedPayoff = max((S - K), 0);
    % % % % %     boundaryMatrixT(1,:) = discountedPayoff;
    % % % % % 
    % % % % %     % Compute first derivative with respect to V
    % % % % %     firstVDerivative = MDerivativeC(u', dv, 1, 1, boundaryMatrixT)';  % Transposed result
    % % % % % 
    % % % % %     % Compute second derivative with respect to V
    % % % % %     secondVDerivative = MSecondDerivativePlusC(u', dv, 1, 1, boundaryMatrixT)';  % Transposed result
    % % % % % 
    % % % % %     % Compute mixed derivative
    % % % % %     boundaryMatrix_mixed = zeros(NS, NV);
    % % % % %     boundaryMatrix_mixed(NS,:) = (1 / ds);  % Neumann condition for mixed derivative
    % % % % %     mixedDerivative = MDerivativeC(firstVDerivative, ds, 0, 1, boundaryMatrix_mixed);
    % % % % % 
    % % % % %     %the following is needed for operator part A3
    % % % % %     IDNV = eye( NV ,'like', VMatrix );
    % % % % % 
    % % % % %     % Compute terms A1 to A5
    % % % % %     A1 = (r - q) * SMatrix * firstSDerivative;
    % % % % %     A2 = 0.5 * (S2Matrix * secondSDerivative) * VMatrix;
    % % % % %     A3 = ((kappa * (theta*IDNV - VMatrix) - lambda * VMatrix) * firstVDerivative')';
    % % % % %     A4 = 0.5 * sigma^2 * (secondVDerivative) * VMatrix;
    % % % % %     A5 = rho * sigma * SMatrix * mixedDerivative * VMatrix;
    % % % % % 
    % % % % %     U = (1-dt*r)*u+dt*(A1+A2+A3+A4+A5);
    % % % % % end


    for t = 1:NT-1

        boundaryMatrix_S = zeros(NS, NV);
        boundaryMatrix_S(NS,:) = (1 / ds);  %Neumann condition at large S

        boundaryMatrix_V = zeros(NS, NV);
        boundaryMatrixT = boundaryMatrix_V';  % Transpose for volatility derivatives
        boundaryMatrixT(NV,:) = S;  % Boundary at max volatility

        % Boundary Condition when V = 0 (volatility = 0)
        discountedPayoff = max((S - K * exp(-r * (Tmax - T(t)))), 0);
        %discountedPayoff = max((S - K), 0);
        boundaryMatrixT(1,:) = discountedPayoff;

        boundaryMatrix_mixed = zeros(NS, NV);
        boundaryMatrix_mixed(NS,:) = (1 / ds);  % Neumann condition for mixed derivative


        % Compute first and second derivatives in S
        firstSDerivative = MDerivativeC(U, ds, 0, 1, boundaryMatrix_S);
        secondSDerivative = MSecondDerivativePlusC(U, ds, 0, 1, boundaryMatrix_S);
        
        % Compute first and second derivatives in V
        firstVDerivative = MDerivativeC(U', dv, 1, 1, boundaryMatrixT)';
        secondVDerivative = MSecondDerivativePlusC(U', dv, 1, 1, boundaryMatrixT)';
    
        % Compute mixed derivative
        mixedDerivative = MDerivativeC(firstVDerivative, ds, 0, 1, boundaryMatrix_mixed);
    
        % Compute A1 to A5 using the derivatives
        A1 = (r - q) * SMatrix * firstSDerivative;
        A2 = 0.5 * (S2Matrix * secondSDerivative) * VMatrix;
        A3 = ((kappa * (theta * eye(NV) - VMatrix) - lambda * VMatrix) * firstVDerivative')';
        A4 = 0.5 * sigma^2 * (secondVDerivative) * VMatrix;
        A5 = rho * sigma * SMatrix * mixedDerivative * VMatrix;
    
        % Combine into matrix A
        A = A1 + A2 + A3 + A4 + A5;
    
        % Crank-Nicolson update
        


        % % % lhs_matrix = (eye(NS, NV) - dt / 2 * A);
        % % % 
        % % % U_vec = U(:);  
        % % % 
        % % % rhs_matrix = (eye(NS, NV) + dt / 2 * A) * U_vec + dt * boundaryMatrix_S(:);
        % % % U = lhs_matrix \ rhs_matrix;

        % % % % Flatten U into a column vector
        % % % U_vec = U(:);  
        % % % 
        % % % % Reshape A into a block matrix
        % % % A_large = kron(eye(NS), A);  % Create block matrix
        % % % 
        % % % % Crank-Nicolson update
        % % % lhs_matrix = (eye(NS*NV) - dt / 2 * A_large);  % LHS matrix
        % % % rhs_matrix = (eye(NS*NV) + dt / 2 * A_large) * U_vec + dt * boundaryMatrix_S(:);  % RHS matrix
        % % % 
        % % % % Solve for U_vec at the next time step
        % % % U_vec = lhs_matrix \ rhs_matrix;
        % % % 
        % % % % Reshape U_vec back to the original NS x NV dimensions
        % % % U = reshape(U_vec, [NS, NV]);



        %%%%%%%% Pedestrian Crank Nicolson
        % % % % % % 
        % % % % % % % Flatten U into a column vector (size: NS*NV x 1)
        % % % % % % U_vec = U(:);
        % % % % % % 
        % % % % % % A_vec = A(:);
        % % % % % % 
        % % % % % % % Perform Crank-Nicolson update without Kronecker product
        % % % % % % lhs_matrix = (eye(NS*NV) - dt / 2 * A_vec);  % Ensure A has size NS*NV x NS*NV
        % % % % % % rhs_matrix = (eye(NS*NV) + dt / 2 * A_vec) * U_vec + dt * boundaryMatrix_S(:);
        % % % % % % 
        % % % % % % % Solve for U_vec at the next time step
        % % % % % % U_vec = lhs_matrix \ rhs_matrix;
        % % % % % % 
        % % % % % % % Reshape U_vec back to the original NS x NV dimensions
        % % % % % % U = reshape(U_vec, [NS, NV]);


        %U is size NS,NV, the vectorised form is NS*NV,1
        U_vec = U(:);

        A_vec = A(:);

        combined_boundary_conditions = boundaryMatrix_S(:) + boundaryMatrixT(:) + boundaryMatrix_mixed(:);
        %combined_boundary_conditions = boundaryMatrix_S(:);

        % LHS and RHS matrices for Crank-Nicolson (NS*NV x NS*NV identity)
        lhs_matrix = ((1-(dt*r/2))*eye(NS*NV) + (dt/2)*A_vec);  % Matrix multiplication directly on A, ensuring it's a vectorized operation
        rhs_matrix = ((1+(dt*r/2))*eye(NS*NV) - (dt/2)*A_vec) * U_vec - dt*combined_boundary_conditions(:);  % Ensure U_vec and boundary conditions are vectors
        
        % Set initial guess for GMRES
        x0 = U_vec;  % The current solution vector
        
        % Set GMRES parameters
        restart = 80;  % Restart after 20 iterations (example value)
        tol = 1e-5;  % Tolerance for convergence
        max_iter = 100;  % Maximum number of iterations
        
        % Solve using GMRES
        warning('off', 'all')
        %U_vec = restarted_gmres(lhs_matrix, rhs_matrix, x0, restart, tol, max_iter);
        U_vec = restarted_gmres(lhs_matrix, rhs_matrix, x0, restart, tol, max_iter);
        warning('on', 'all')
        % Reshape U_vec back to the original NS x NV dimensions
        U = reshape(U_vec, [NS, NV]);
    end    

end

function R = SVDRank(U, epsilon)

    %builtin_rank = rank(U,epsilon);

    [uc,sc,vs]=svd(U);

    sigSum = 0;

    %the look should go to the minimum of the two dimensions
    %so a step is missing and it should be added when moving
    %this to a separate function
    for i = 1: size(sc,2)
        sigSum = sigSum + power(sc(i,i),2);
    end

    threshold = power(epsilon,2)*sigSum;

    rank_home_made = 0;

    errorSquare = 0;

    for i=size(sc,2):-1:1
        errorSquare = errorSquare + power(sc(i,i),2);
        if errorSquare > threshold
            rank_home_made = i + 1;
            break;
        end
    end

    R=rank_home_made;
end

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


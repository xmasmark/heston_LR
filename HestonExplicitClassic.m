function U = HestonExplicitClassic(params,K,r,q,S,V,T)

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
    for t = 1:NT-1
        % Initialize u for this time step
        u = U;

        % Set up boundary matrices and derivatives in S
        boundaryMatrix_S = zeros(NS, NV);
        boundaryMatrix_S(NS,:) = (1 / ds);  %Neumann condition at large S

        firstSDerivative = MDerivativeC(u, ds, 0, 1, boundaryMatrix_S);
        secondSDerivative = MSecondDerivativePlusC(u, ds, 0, 1, boundaryMatrix_S);

        % Set up boundary matrices and derivatives in V
        boundaryMatrix_V = zeros(NS, NV);
        boundaryMatrixT = boundaryMatrix_V';  % Transpose for volatility derivatives
        boundaryMatrixT(NV,:) = S;  % Boundary at max volatility

        % Boundary Condition when V = 0 (volatility = 0)
        discountedPayoff = max((S - K * exp(-r * (Tmax - T(t)))), 0);
        %discountedPayoff = max((S - K), 0);
        boundaryMatrixT(1,:) = discountedPayoff;

        % Compute first derivative with respect to V
        firstVDerivative = MDerivativeC(u', dv, 1, 1, boundaryMatrixT)';  % Transposed result

        % Compute second derivative with respect to V
        secondVDerivative = MSecondDerivativePlusC(u', dv, 1, 1, boundaryMatrixT)';  % Transposed result

        % Compute mixed derivative
        boundaryMatrix_mixed = zeros(NS, NV);
        boundaryMatrix_mixed(NS,:) = (1 / ds);  % Neumann condition for mixed derivative
        mixedDerivative = MDerivativeC(firstVDerivative, ds, 0, 1, boundaryMatrix_mixed);

        %the following is needed for operator part A3
        IDNV = eye( NV ,'like', VMatrix );

        % Compute terms A1 to A5
        A1 = (r - q) * SMatrix * firstSDerivative;
        A2 = 0.5 * (S2Matrix * secondSDerivative) * VMatrix;
        A3 = ((kappa * (theta*IDNV - VMatrix) - lambda * VMatrix) * firstVDerivative')';
        A4 = 0.5 * sigma^2 * (secondVDerivative) * VMatrix;
        A5 = rho * sigma * SMatrix * mixedDerivative * VMatrix;

        U = (1-dt*r)*u+dt*(A1+A2+A3+A4+A5);
    end

    % % % % % % % %attempting GMRES
    % % % % % % % for t = 1:NT-1
    % % % % % % %     % Initialize u for this time step
    % % % % % % %     u = U;
    % % % % % % % 
    % % % % % % %     % Flatten u into a vector for the Crank-Nicolson scheme
    % % % % % % %     u_vec = reshape(u, [], 1);  % Flatten U into a 1D vector (NS * NV x 1)
    % % % % % % % 
    % % % % % % %     % Set up boundary matrices and derivatives (same as before)
    % % % % % % %     boundaryMatrix_S = zeros(NS, NV);
    % % % % % % %     boundaryMatrix_S(NS,:) = (1 / ds);  % Neumann condition at large S
    % % % % % % % 
    % % % % % % %     % Compute first and second derivatives (as before)
    % % % % % % %     firstSDerivative_n = MDerivativeC(u, ds, 0, 1, boundaryMatrix_S);
    % % % % % % %     secondSDerivative_n = MSecondDerivativePlusC(u, ds, 0, 1, boundaryMatrix_S);
    % % % % % % % 
    % % % % % % %     % Compute derivatives with respect to V (as before)
    % % % % % % %     boundaryMatrix_V = zeros(NS, NV);
    % % % % % % %     boundaryMatrixT = boundaryMatrix_V';
    % % % % % % %     boundaryMatrixT(NV,:) = S;
    % % % % % % %     discountedPayoff = max((S - K * exp(-r * (Tmax - T(t)))), 0);
    % % % % % % %     boundaryMatrixT(1,:) = discountedPayoff;
    % % % % % % % 
    % % % % % % %     firstVDerivative_n = MDerivativeC(u', dv, 1, 1, boundaryMatrixT)';
    % % % % % % %     secondVDerivative_n = MSecondDerivativePlusC(u', dv, 1, 1, boundaryMatrixT)';
    % % % % % % % 
    % % % % % % %     % Compute mixed derivative (as before)
    % % % % % % %     boundaryMatrix_mixed = zeros(NS, NV);
    % % % % % % %     boundaryMatrix_mixed(NS,:) = (1 / ds);
    % % % % % % %     mixedDerivative_n = MDerivativeC(firstVDerivative_n, ds, 0, 1, boundaryMatrix_mixed);
    % % % % % % % 
    % % % % % % %     % Compute terms A1 to A5
    % % % % % % %     A1_n = (r - q) * SMatrix * firstSDerivative_n;
    % % % % % % %     A2_n = 0.5 * (S2Matrix * secondSDerivative_n) * VMatrix;
    % % % % % % %     A3_n = ((kappa * (theta * eye(NV, 'like', VMatrix) - VMatrix) - lambda * VMatrix) * firstVDerivative_n')';
    % % % % % % %     A4_n = 0.5 * sigma^2 * (secondVDerivative_n) * VMatrix;
    % % % % % % %     A5_n = rho * sigma * SMatrix * mixedDerivative_n * VMatrix;
    % % % % % % % 
    % % % % % % %     % Right-hand side (explicit terms)
    % % % % % % %     RHS = (1 - dt * r) * u + (dt / 2) * (A1_n + A2_n + A3_n + A4_n + A5_n);
    % % % % % % % 
    % % % % % % %     % Flatten RHS to vector form
    % % % % % % %     RHS_vec = reshape(RHS, [], 1);
    % % % % % % % 
    % % % % % % %     % Now construct the matrix A = I - dt/2 * A_new for GMRES
    % % % % % % %     A_new = A1_n + A2_n + A3_n + A4_n + A5_n;  % Update A_new for implicit part
    % % % % % % %     A_new_flat = reshape(A_new, [], 1);  % Flatten A_new for matrix operations
    % % % % % % % 
    % % % % % % %     LHS = speye(NS * NV) - (dt / 2) * A_new_flat;  % Left-hand side matrix
    % % % % % % % 
    % % % % % % %     % Use GMRES to solve the linear system (LHS * U_new_vec = RHS_vec)
    % % % % % % %     % Initial guess for GMRES (could be the current u_vec or zeros(NS*NV, 1))
    % % % % % % %     x0 = u_vec;  % Initial guess for GMRES
    % % % % % % % 
    % % % % % % %     % Parameters for GMRES
    % % % % % % %     restart = 50;  % Number of iterations before restarting GMRES
    % % % % % % %     tol = 1e-6;    % Tolerance for GMRES
    % % % % % % %     max_iter = 100; % Maximum number of iterations
    % % % % % % % 
    % % % % % % %     % Solve using GMRES
    % % % % % % %     U_new_vec = restarted_gmres(LHS, RHS_vec, x0, restart, tol, max_iter);
    % % % % % % % 
    % % % % % % %     % Reshape U_new_vec back to 2D form
    % % % % % % %     U_new = reshape(U_new_vec, [NS, NV]);
    % % % % % % % 
    % % % % % % %     % Update U for the next iteration
    % % % % % % %     U = U_new;
    % % % % % % % end





    % % % % % for t = 1:NT-1
    % % % % %     % Initialize u for this time step
    % % % % %     u = U;
    % % % % % 
    % % % % %     % Flatten u into a vector for the Crank-Nicolson scheme
    % % % % %     u_vec = reshape(u, [], 1);  % Flatten U into a 1D vector (NS * NV x 1)
    % % % % % 
    % % % % %     % Set up boundary matrices and derivatives in S
    % % % % %     boundaryMatrix_S = zeros(NS, NV);
    % % % % %     boundaryMatrix_S(NS,:) = (1 / ds);  % Neumann condition at large S
    % % % % % 
    % % % % %     % Compute first and second derivatives in S (flattened approach)
    % % % % %     firstSDerivative_n = MDerivativeC(u, ds, 0, 1, boundaryMatrix_S);
    % % % % %     secondSDerivative_n = MSecondDerivativePlusC(u, ds, 0, 1, boundaryMatrix_S);
    % % % % % 
    % % % % %     % Set up boundary matrices and derivatives in V
    % % % % %     boundaryMatrix_V = zeros(NS, NV);
    % % % % %     boundaryMatrixT = boundaryMatrix_V';  % Transpose for volatility derivatives
    % % % % %     boundaryMatrixT(NV,:) = S;  % Boundary at max volatility
    % % % % % 
    % % % % %     % Boundary Condition when V = 0 (volatility = 0)
    % % % % %     discountedPayoff = max((S - K * exp(-r * (Tmax - T(t)))), 0);
    % % % % %     boundaryMatrixT(1,:) = discountedPayoff;
    % % % % % 
    % % % % %     % Compute first and second derivatives with respect to V (flattened approach)
    % % % % %     firstVDerivative_n = MDerivativeC(u', dv, 1, 1, boundaryMatrixT)';  % Transposed result
    % % % % %     secondVDerivative_n = MSecondDerivativePlusC(u', dv, 1, 1, boundaryMatrixT)';  % Transposed result
    % % % % % 
    % % % % %     % Compute mixed derivative
    % % % % %     boundaryMatrix_mixed = zeros(NS, NV);
    % % % % %     boundaryMatrix_mixed(NS,:) = (1 / ds);  % Neumann condition for mixed derivative
    % % % % %     mixedDerivative_n = MDerivativeC(firstVDerivative_n, ds, 0, 1, boundaryMatrix_mixed);
    % % % % % 
    % % % % %     % Needed for operator part A3
    % % % % %     IDNV = eye(NV, 'like', VMatrix);
    % % % % % 
    % % % % %     % Compute terms A1 to A5 at the current time step (explicit part)
    % % % % %     A1_n = (r - q) * SMatrix * firstSDerivative_n;
    % % % % %     A2_n = 0.5 * (S2Matrix * secondSDerivative_n) * VMatrix;
    % % % % %     A3_n = ((kappa * (theta * IDNV - VMatrix) - lambda * VMatrix) * firstVDerivative_n')';
    % % % % %     A4_n = 0.5 * sigma^2 * (secondVDerivative_n) * VMatrix;
    % % % % %     A5_n = rho * sigma * SMatrix * mixedDerivative_n * VMatrix;
    % % % % % 
    % % % % %     % Right-hand side (explicit terms)
    % % % % %     RHS = (1 - dt * r) * u + (dt / 2) * (A1_n + A2_n + A3_n + A4_n + A5_n);
    % % % % % 
    % % % % %     % Flatten RHS to vector form
    % % % % %     RHS_vec = reshape(RHS, [], 1);  % Flatten RHS to 1D vector
    % % % % % 
    % % % % %     % Now compute the implicit part (future time step)
    % % % % %     % The operator at time step n+1 will involve U at the next time step, so we need to solve
    % % % % %     % U_new = RHS + dt/2 * A_new
    % % % % %     % Placeholder for the operators at the next step
    % % % % %     A1_new = (r - q) * SMatrix * firstSDerivative_n;  % Placeholder, needs update
    % % % % %     A2_new = 0.5 * (S2Matrix * secondSDerivative_n) * VMatrix;  % Placeholder
    % % % % %     A3_new = ((kappa * (theta * IDNV - VMatrix) - lambda * VMatrix) * firstVDerivative_n')';  % Placeholder
    % % % % %     A4_new = 0.5 * sigma^2 * (secondVDerivative_n) * VMatrix;  % Placeholder
    % % % % %     A5_new = rho * sigma * SMatrix * mixedDerivative_n * VMatrix;  % Placeholder
    % % % % % 
    % % % % %     % Combine for Crank-Nicolson update
    % % % % %     A_new = A1_new + A2_new + A3_new + A4_new + A5_new;
    % % % % % 
    % % % % %     % Flatten A_new to handle in vectorized form (assuming you compute A_new as a sparse matrix)
    % % % % %     A_new_flat = reshape(A_new, [], 1);  % Flatten the new operator
    % % % % % 
    % % % % %     % Construct the linear system for the Crank-Nicolson update
    % % % % %     LHS = speye(NS * NV) - (dt / 2) * A_new_flat;  % Use sparse identity matrix
    % % % % % 
    % % % % %     % Solve for U at the next time step using a linear solver
    % % % % %     U_new_vec = LHS \ RHS_vec;  % Solve the linear system in vectorized form
    % % % % % 
    % % % % %     % Reshape U_new back to a matrix after solving
    % % % % %     U_new = reshape(U_new_vec, [NS, NV]);
    % % % % % 
    % % % % %     % Update U for the next iteration
    % % % % %     U = U_new;
    % % % % % end
end


%writematrix(matrix_rank, 'G:\My Drive\BO\PhD\Heston\RankAnalysis\matrix_rank.csv')
%writematrix(matrix_rank_home_made, 'G:\My Drive\BO\PhD\Heston\RankAnalysis\matrix_rank_home_made.csv')

%Task: rewrite the code for 1 step only (t=0->t=1)
%In the new code I don't use 2 dimensional U and u: it has to work with
%2 vectors that are 1 dimensional (1 row and 1 column).
%the start is U that is 2 vectors, 1 depending on S and one depending on V.
%Like it was done for the boundaries. The stuff that is currently written
%in the two dimensional case, needs to work on vectors rather than 2
%dimensional matrix. One step implementation will show the linear
%complexity leading to the same result of the quadratic complexity (current
%implementation).

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

%for now it is better to pre-compute the boundaries in the code above and
%to drive this method through the precomputed boundaries
function D = MDerivative(matrix, delta, boundaryConditionLeft, boundaryConditionRight, boundaryMatrix, direction)
    %boundaryCondition -- Dirichlet = 0, Neumann = 1
    %boundaryMatrix are the actual correction factors to apply
    %to insert in the correction vector
    %direction is an important parameter: 0 and 1. 
    %In the case of zero it derivates by rows
    %In the case of one it derivates by columns
    %this is because the heston matrix is composed by two variables: S and V
    %... mind you, I could just transpose before calculating the derivative
    %and that would save expensive lines of code.

    % Calculate the size of the input matrix L
    [records, fields] = size(matrix);
   
    direction = 0;

    n = records;
    h = delta;

    %boundaryVector = zeros(n,1);

    % Create the "Laplacian" matrix
    temporary = zeros(n, n);
    temporary = temporary - diag(ones(n-1, 1), -1);
    temporary = temporary + diag(ones(n-1, 1), 1);
    
    % % Handle boundary conditions
    if(boundaryConditionLeft == 0)
        temporary(1,1)=0;
    elseif(boundaryConditionLeft == 1)
        temporary(1,1)=-1;
    end

    if(boundaryConditionRight == 0)
        temporary(n,n)=0;
    elseif(boundaryConditionRight == 1)
        temporary(n,n)=1;
    end

    % Scale by 1/(2*h) to account for the step size
    temporary = temporary / (2 * h);

    % Compute the derivative
    D = (temporary * matrix) + boundaryMatrix;

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


function D = MSecondDerivative(L, delta)
    % Calculate the size of the input matrix L
    [records, fields] = size(L);

    n = records;
    h = delta;
    matrix = L;

    % Create the Laplacian matrix for second derivative
    laplacian = zeros(n, n);
    laplacian = laplacian + diag(ones(n-1, 1), -1);  % sub-diagonal
    laplacian = laplacian - 2 * diag(ones(n, 1), 0); % main diagonal
    laplacian = laplacian + diag(ones(n-1, 1), 1);   % super-diagonal
    
    % Handle boundary conditions (simple forward and backward second difference)
    laplacian(1, 1) = 1;
    laplacian(1, 2) = -2;
    laplacian(1, 3) = 1;
    laplacian(n, n-2) = 1;
    laplacian(n, n-1) = -2;
    laplacian(n, n) = 1;

    % Scale by 1/h^2 to account for the step size
    laplacian = laplacian / (h^2);

    % Compute the second derivative
    D = laplacian * matrix;    
end


function D = MSecondDerivativePlus(L, delta, boundaryConditionLeft, boundaryConditionRight, boundaryMatrix)
    % Calculate the size of the input matrix L
    [records, fields] = size(L);

    n = records;
    h = delta;
    matrix = L;

    % Create the Laplacian matrix for second derivative
    laplacian = zeros(n, n);
    laplacian = laplacian + diag(ones(n-1, 1), -1);  % sub-diagonal
    laplacian = laplacian - 2 * diag(ones(n, 1), 0); % main diagonal
    laplacian = laplacian + diag(ones(n-1, 1), 1);   % super-diagonal
    
    % Handle boundary conditions (simple forward and backward second difference)
    laplacian(1, 1) = boundaryConditionLeft;
    % laplacian(1, 2) = -2;
    % laplacian(1, 3) = 1;
    % laplacian(n, n-2) = 1;
    % laplacian(n, n-1) = -2;
    laplacian(n, n) = boundaryConditionRight;

    % Scale by 1/h^2 to account for the step size
    laplacian = laplacian / (h^2);

    % Compute the second derivative
    D = laplacian * matrix + boundaryMatrix;    
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


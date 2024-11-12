function U = HestonExplicitClassicCNXYRC01(params,K,r,q,S,V,T,mode)

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

        %I leave the following two statements
        x = [X];
        y = [Y];

        boundaryMatrix_S = zeros(NS, NV);
        boundaryMatrix_S(NS,:) = (1/ds);  %Neumann condition at large S

        bm_S=boundaryMatrix_S;

        b1x = zeros(NS,1);
        b1x(NS)=(1/ds);
        b1y = ones(NV,1);


        boundaryMatrix_S2 = zeros(NS, NV);
        boundaryMatrix_S2(NS,:) = (1/ds);  %Neumann condition at large S

        bm_S2=boundaryMatrix_S2;

        b2x = zeros(NS,1);
        b2x(NS)=(1/ds);
        b2y = ones(NV,1);
        
        boundaryMatrix_V = zeros(NS, NV);
        boundaryMatrixT = boundaryMatrix_V';  % Transpose for volatility derivatives
        boundaryMatrixT(NV,:) = S;  % Boundary at max volatility

        % Boundary Condition when V = 0 (volatility = 0)
        discountedPayoff = max((S - K * exp(-r * (Tmax - T(t)))), 0);
        %discountedPayoff = max((S - K), 0);
        boundaryMatrixT(1,:) = discountedPayoff;

        bm_V1 = boundaryMatrixT;

        b3x = [discountedPayoff', S'];%In this case rank 2 because of the shape of the condition
        b3y = zeros(NV,2); %2 columns because this is rank 2
        b3y(1,1)=1;
        b3y(NV,2)=1;

        boundaryMatrix_mixed = zeros(NS, NV);
        boundaryMatrix_mixed(NS,:) = (1 / ds);  % Neumann condition for mixed derivative

        bm_SV = boundaryMatrix_mixed;

        d1sM = MDerivativeVM(U, ds, 0, 1);
        d2sM = MSecondDerivativePlusCVM(U, ds, 0, 1);
        
        b5x = zeros(NS,1);
        b5x(NS)=(1/ds);
        b5y = ones(NV,1);

        % Compute first and second derivatives in V
    
        d1vM = MDerivativeVM(U', dv, 1, 1);
        d2vM = MSecondDerivativePlusCVM(U', dv, 1, 1);

        IS = eye(NS);
        IV = eye(NV);

        A1KX = (r-q)*SMatrix*d1sM;
        A1KY = IV;
        A2KX = 0.5*S2Matrix*d2sM;
        A2KY = VMatrix;
        A3KX = IS;
        A3KY = (kappa * (theta * eye(NV) - VMatrix) - lambda * VMatrix)* d1vM;
        A4KX = IS;
        A4KY = (0.5 * sigma^2)*VMatrix * d2vM;
        A5KX = SMatrix*d1sM;
        A5KY = rho * sigma * (VMatrix*d1vM);

        X1 = A1KX*x;
        Y1 = A1KY*y;
        X2 = A2KX*x;
        Y2 = A2KY*y;
        X3 = A3KX*x;
        Y3 = A3KY*y;
        X4 = A4KX*x;
        Y4 = A4KY*y;
        X5 = A5KX*x;
        Y5 = A5KY*y;

        b1x = (r-q) * SMatrix * b1x;
        b2x = 0.5 * S2Matrix * b2x;
        b2y = VMatrix'*b2y;
        b3y = (kappa * (theta * eye(NV) - VMatrix) - lambda * VMatrix)' * b3y;
        b4x = zeros(NS,1);
        b4y = zeros(NV,1);
        
        b4y = 0.5 * sigma^2*VMatrix'*b4y;
        b5x = rho * sigma * SMatrix * b5x;
        b5y = VMatrix'*b5y;

        % AKX = [A1KX, A2KX, A3KX, A4KX, A5KX];
        % AKY = [A1KY, A2KY, A3KY, A4KY, A5KY];
        AX = [X1, X2, X3, X4, X5];
        AY = [Y1, Y2, Y3, Y4, Y5];


        % % % %     U_vec = (1-dt*r)*U_vec+dt*(Av + Bv);


        % % the following applies in the case of Euler
        % % dtsrt = sqrt(dt);
        % % 
        % % XC = dtsrt*[X1,X2,X3,X4,X5,b1x,b2x,b3x,b4x,b5x];
        % % YC = dtsrt*[Y1,Y2,Y3,Y4,Y5,b1y,b2y,b3y,b4y,b5y];
        % % 
        % % cs = sqrt((1-dt*r));
        % % 
        % % newX = [cs*x XC];
        % % newY = [cs*y YC];

        % Set GMRES parameters
        restart = 80;  % Restart after 20 iterations (example value)
        tol = 1e-5;  % Tolerance for convergence
        max_iter = 100;  % Maximum number of iterations

        % x0=zeros(NS,size(AKX,2));
        % y0=zeros(NV,size(AKY,2));
        % x0=zeros(NS,size(AX,2));
        % y0=zeros(NV,size(AY,2));
        % x0=zeros(NS,NS);
        % y0=zeros(NV,NV);

        % XC = [X1,X2,X3,X4,X5,b1x,b2x,b3x,b4x,b5x];
        % YC = [Y1,Y2,Y3,Y4,Y5,b1y,b2y,b3y,b4y,b5y];

        [X_new, Y_new] = GMRES_XYv01(AX, AY, x, y, restart, tol, max_iter);
        %[X_new, Y_new] = GMRES_XYv01(AX, AY, XC, YC, restart, tol, max_iter);

        [X,Y]=CompressData(X_new, Y_new,epsilon);


        %fprintf('rank of XC is %d\n', rank(XC));

        %[X,Y]=CompressData(newX,newY,epsilon);

        

        %%%Comment:
        %%%before this step I need to plug in an efficient GMRES
        %%%implementation

        % % % % lhs_matrix = ((1+(dt*r/2))*eye(NS*NV) - (dt/2)*AK);  
        % % % % rhs_vector = ((1-(dt*r/2))*U_vec + (dt/2)*Av) + dt*Bv;
        % AK in low rank format: AKx*Aky'
        % AV in low rank format: 


        % % % % %%%OLD CODE STARTS
        % % % % AK= A1k+A2k+A3k+A4k+A5k;
        % % % % 
        % % % % A1v = A1k * U_vec;
        % % % % A2v = A2k * U_vec;
        % % % % A3v = A3k * U_vec;
        % % % % A4v = A4k * U_vec;
        % % % % A5v = A5k * U_vec;
        % % % % 
        % % % % Av = A1v+A2v+A3v+A4v+A5v;
        % % % % 
        % % % % B1v = (r-q) * SMatrix * d1SB;
        % % % % B2v = 0.5 * S2Matrix * d2SB * VMatrix;
        % % % % B3v = ((kappa * (theta * eye(NV) - VMatrix) - lambda * VMatrix) * d1VB)';
        % % % % B4v = 0.5 * sigma^2 * d2VB' * VMatrix;
        % % % % B5v = rho * sigma * SMatrix * dSVB * VMatrix;
        % % % % 
        % % % % Bv = B1v(:)+B2v(:)+B3v(:)+B4v(:)+B5v(:);
        % % % % 
        % % % % lhs_matrix = ((1+(dt*r/2))*eye(NS*NV) - (dt/2)*AK);  
        % % % % rhs_vector = ((1-(dt*r/2))*U_vec + (dt/2)*Av) + dt*Bv;
        % % % % 
        % % % % % Solve for U_vec at the next time step
        % % % % if mode == 0
        % % % %     U_vec = (1-dt*r)*U_vec+dt*(Av + Bv);
        % % % % end
        % % % % if mode == 1
        % % % %     U_vec = lhs_matrix \ rhs_vector;
        % % % % end
        % % % % if mode == 2
        % % % %     % Set initial guess for GMRES
        % % % %     x0 = U_vec;  % The current solution vector
        % % % % 
        % % % %     % Set GMRES parameters
        % % % %     restart = 80;  % Restart after 20 iterations (example value)
        % % % %     tol = 1e-5;  % Tolerance for convergence
        % % % %     max_iter = 100;  % Maximum number of iterations
        % % % % 
        % % % %     % Solve using GMRES
        % % % %     warning('off', 'all')
        % % % %     U_vec = restarted_gmres(lhs_matrix, rhs_vector, x0, restart, tol, max_iter);
        % % % %     warning('on', 'all')
        % % % % end
        % % % % %%%OLD CODE ENDS
        % % % % 


        %U = (1-dt*r)*U+dt*(A1+A2+A3+A4+A5);

        %this is in matrix form and works beautifully
        %U = (1-dt*r)*U+dt*(A + B);

        %this is in vector form and works beautifully
        %U_vec_Euler = (1-dt*r)*U_vec+dt*(Av + Bv);

        % % % % % % % if mod(t,100)==0
        % % % % % % %     condition = cond(lhs_matrix);
        % % % % % % %     fprintf('condition at time step %d is %d\n', t, condition);
        % % % % % % %     diff = reshape(U_vec_Euler, [NS, NV]) - reshape(U_vec_Crank_Nicolson, [NS, NV]);
        % % % % % % %     fprintf('The norm of the difference at time step %d is %d\n', t, norm(diff));
        % % % % % % % end

        % Reshape U_vec back to the original NS x NV dimensions
        %U = reshape(U_vec_Euler, [NS, NV]);


        %%%%%%%U = reshape(U_vec, [NS, NV]);

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


function U = HestonExplicitClassicCNXYCOMP(params,K,r,q,S,V,T)

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

    testDiff = U-X*Y';

    for t = 1:NT-1

        %here is the low rank part
        tol = 1e-3;  % Tolerance for convergence and compression
        [x,y]=CompressData(X,Y,tol);
        [AX,AY] = HestonMatVec(x,y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);
        [BX,BY] = HestonMatVecBoundaries(NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, K, Tmax, t, T);

        %half Euler step
        FX = [(1+r*dt/2)*x, (-dt/2)*AX, BX]; 
        FY = [           y,         AY, BY];

        %Right hand side vector components
        [BXc,BYc]=CompressData(FX, FY, epsilon);
        %the LHS are operators, not a matrix
        % Set GMRES parameters
        restart = 3;  % Restart after N iterations
        max_iter = 10;  % Maximum number of iterations

        %%x and y, old values, the initial guesses
        %[X, Y] = GMRES_XYv01(x, y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, K, Tmax, t, T, BXc, BYc, x, y, restart, tol, max_iter);

        %here is the vectorised part

        boundaryMatrix_S = zeros(NS, NV);
        boundaryMatrix_S(NS,:) = (1/ds);  %Neumann condition at large S
        bm_S=boundaryMatrix_S;

        boundaryMatrix_S2 = zeros(NS, NV);
        boundaryMatrix_S2(NS,:) = (1/ds);  %Neumann condition at large S
        bm_S2=boundaryMatrix_S2;

        boundaryMatrix_V = zeros(NS, NV);
        boundaryMatrixT = boundaryMatrix_V';  % Transpose for volatility derivatives
        boundaryMatrixT(NV,:) = S;  % Boundary at max volatility

        % Boundary Condition when V = 0 (volatility = 0)
        discountedPayoff = max((S - K * exp(-r * (Tmax - T(t)))), 0);
        %discountedPayoff = max((S - K), 0);
        boundaryMatrixT(1,:) = discountedPayoff;

        bm_V1 = boundaryMatrixT;

        boundaryMatrix_mixed = zeros(NS, NV);
        boundaryMatrix_mixed(NS,:) = (1 / ds);  % Neumann condition for mixed derivative

        bm_SV = boundaryMatrix_mixed;

        d1sM = MDerivativeVM(NS, ds, 0, 1);
        d2sM = MSecondDerivativePlusCVM(NS, ds, 0, 1);
        
        d1SB = bm_S;
        d2SB = bm_S2;
        d1VB = bm_V1;
        d2VB = zeros(NV,NS);
        dSVB = bm_SV;

        % Compute first and second derivatives in V
    
        d1vM = MDerivativeVM(NV, dv, 1, 1);
        d2vM = MSecondDerivativePlusCVM(NV, dv, 1, 1);

        IS = eye(NS);
        IV = eye(NV);

        %U is size NS,NV, the vectorised form is NS*NV,1
        U_vec = U(:);

        A1k = (r-q)*kron(IV,SMatrix*d1sM);
        A2k = 0.5 * kron(VMatrix,S2Matrix*d2sM);
        A3k = (kron((kappa * (theta * eye(NV) - VMatrix) - lambda * VMatrix)* d1vM,IS));
        A4k = 0.5 * sigma^2 * kron(VMatrix * d2vM, IS);
        A5k = rho * sigma * kron(VMatrix*d1vM, SMatrix*d1sM);

        AK= A1k+A2k+A3k+A4k+A5k;

        A1v = A1k * U_vec;
        A2v = A2k * U_vec;
        A3v = A3k * U_vec;
        A4v = A4k * U_vec;
        A5v = A5k * U_vec;

        Av = A1v+A2v+A3v+A4v+A5v;

        B1v = (r-q) * SMatrix * d1SB;
        B2v = 0.5 * S2Matrix * d2SB * VMatrix;
        B3v = ((kappa * (theta * eye(NV) - VMatrix) - lambda * VMatrix) * d1VB)';
        B4v = 0.5 * sigma^2 * d2VB' * VMatrix;
        B5v = rho * sigma * SMatrix * dSVB * VMatrix;

        Bv = B1v(:)+B2v(:)+B3v(:)+B4v(:)+B5v(:);

        lhs_matrix = ((1+(dt*r/2))*eye(NS*NV) - (dt/2)*AK);  
        rhs_vector = ((1-(dt*r/2))*U_vec + (dt/2)*Av) + dt*Bv;

        % Solve for U_vec at the next time step
        % Set initial guess for GMRES
        x0v = U_vec;  % The current solution vector

        % Set GMRES parameters
        restart = 3;  % Restart after 20 iterations (example value)
        tol = 1e-3;  % Tolerance for convergence
        max_iter = 10;  % Maximum number of iterations
            
        % Solve using GMRES
        warning('off', 'all')
        %U_vec = restarted_gmres(lhs_matrix, rhs_vector, x0, restart, tol, max_iter);
        [X_new, Y_new] = GMRES_XY_COMP(lhs_matrix, rhs_vector, x0v, x, y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, K, Tmax, t, T, BXc, BYc, x, y, restart, tol, max_iter);
        warning('on', 'all')
        
        U = reshape(U_vec, [NS, NV]);
        reconstructedU = X*Y';

        diff = U-reconstructedU;

        fprintf('Norm of diff %g\n', norm(diff, 'fro'));

    end    
    U=X*Y';
end






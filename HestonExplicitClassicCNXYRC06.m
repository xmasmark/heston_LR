function [U, SNU, VNU] = HestonExplicitClassicCNXYRC06(params,K,r,q,S,V,T, NT, iterations, restart)

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
    nnT = length(T);
    Smin = S(1);  Smax = S(NS);
    Vmin = V(1);  Vmax = V(NV);
    Tmin = T(1);  Tmax = T(nnT);
    

    % dS = diff(S_cheb);  % Length NS
    % dV = diff(V_cheb);  % Length NV


    % Increment for Stock Price, Volatility, and Maturity
    ds = (Smax-Smin)/(NS-1);
    dv = (Vmax-Vmin)/(NV-1);
    dt = (Tmax-Tmin)/(nnT-1);
    
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

    % X=zeros(NS,1);
    % for s=1:NS
    %     X(s)=max(S(s)-K,0);
    % end
    % 
    % Y=zeros(NV,1);
    % for v=1:NV
	%     Y(v)=1;
    % end

    NS1 = NS + 1;
    NV1 = NV + 1;


    dsv = ChebyshevGrid(NS, Smin, Smax);   % length = NS+1
    dvv = ChebyshevGrid(NV, Vmin, Vmax);   % length = NV+1

    SNU = dsv;
    VNU = dvv;

    deltaS = diff(dsv);  % length = NS
    deltaV = diff(dvv);  % length = NV
    
    X=zeros(NS1,1);
    for s=1:NS1
        X(s)=max(dsv(s)-K,0);
    end

    Y=zeros(NV1,1);
    for v=1:NV1
	    Y(v)=1;
    end
    
    %the following calculations are always the same at each cycle.
    %therefore they have been taken outside the loop and cached as
    %parameters

    % dsv = ChebyshevGrid(NS,Smin,Smax);
    % dvv = ChebyshevGrid(NV,Vmin,Vmax);
    % 
    % deltaS = diff(dsv);
    % deltaV = diff(dvv);
    % 
    % % SMatrix = diag(S);
    % % S2Matrix = diag(S.^2);
    % % VMatrix = diag(V);
    % 
    % SMatrix = diag(dsv);
    % S2Matrix = diag(dsv.^2);
    % VMatrix = diag(dvv);
    % 
    % 
    % d1sM = MDerivativeVMNU(NS, deltaS, 0, 1);
    % d2sM = MSecondDerivativePlusCVMNU(NS, deltaS, 0, 1);
    % 
    % d1vM = MDerivativeVMNU(NV, deltaV, 1, 1);
    % d2vM = MSecondDerivativePlusCVMNU(NV, deltaV, 1, 1);


    SMatrix  = diag(dsv);
    S2Matrix = diag(dsv.^2);
    VMatrix  = diag(dvv);

    d1sM = MDerivativeVMNU(NS1, deltaS, 0, 1);         % Use NS1 = NS + 1
    d2sM = MSecondDerivativePlusCVMNU(NS1, deltaS, 0, 1);

    d1vM = MDerivativeVMNU(NV1, deltaV, 1, 1);         % Use NV1 = NV + 1
    d2vM = MSecondDerivativePlusCVMNU(NV1, deltaV, 1, 1);
    

    % IS = eye(NS);
    % IV = eye(NV);
    
    IS = eye(NS1);
    IV = eye(NV1);

    A1KX = (r-q)*SMatrix*d1sM;
    A1KY = IV;
    A2KX = 0.5*S2Matrix*d2sM;
    A2KY = VMatrix;
    A3KX = IS;
    % A3KY = (kappa * (theta * eye(NV) - VMatrix) - lambda * VMatrix)* d1vM;
    A3KY = (kappa * (theta * eye(NV1) - VMatrix) - lambda * VMatrix)* d1vM;    
    A4KX = IS;
    A4KY = (0.5 * sigma^2)*VMatrix * d2vM;
    A5KX = SMatrix*d1sM;
    A5KY = rho * sigma * (VMatrix*d1vM);

    % SMatrix = diag(S);
    % S2Matrix = diag(S.^2);
    % VMatrix = diag(V);

    % b1x = zeros(NS,1);
    % b1x(NS)=(1/ds);
    % b1y = ones(NV,1);
    % 
    % b2x = zeros(NS,1);
    % b2x(NS)=(1/ds);
    % b2y = ones(NV,1);
    % 
    % b3y = zeros(NV,2); %2 columns because this is rank 2
    % b3y(1,1)=1;
    % b3y(NV,2)=1;

    % b5x = zeros(NS1,1);
    % b5x(NS1)=(1/ds);
    % b5y = ones(NV,1);

    % b1x = zeros(NS1,1);
    % b1x(NS)=(1/ds);
    % b1y = ones(NV1,1);
    % 
    % b2x = zeros(NS1,1);
    % b2x(NS)=(1/ds);
    % b2y = ones(NV1,1);
    % 
    % b3y = zeros(NV1,2); %2 columns because this is rank 2
    % b3y(1,1)=1;
    % b3y(NV1,2)=1;
    % 
    % b5x = zeros(NS1,1);
    % b5x(NS1)=(1/ds);
    % b5y = ones(NV1,1);

    b1x = zeros(NS1,1);
    b1x(end) = 1 / deltaS(end);  % spacing at the last interval
    b1y = ones(NV1,1);

    b2x = zeros(NS1,1);
    b2x(end) = 1 / deltaS(end);
    b2y = ones(NV1,1);

    b3y = zeros(NV1,2);  % 2 columns → rank-2 correction
    b3y(1,1) = 1;        % bottom boundary in V
    b3y(end,2) = 1;      % top boundary in V

    b5x = zeros(NS1,1);
    b5x(end) = 1 / deltaS(end);  % again, last spacing
    b5y = ones(NV1,1);

    b1x = (r-q) * SMatrix * b1x;
    b2x = 0.5 * S2Matrix * b2x;
    b2y = VMatrix'*b2y;
    b3y = (kappa * (theta * eye(NV1) - VMatrix) - lambda * VMatrix)' * b3y;
    b4x = zeros(NS1,1);
    b4y = zeros(NV1,1);
    
    b4y = 0.5 * sigma^2*VMatrix'*b4y;
    b5x = rho * sigma * SMatrix * b5x;
    b5y = VMatrix'*b5y;

    [A,B] = HestonModelOperatorLean(A1KX,A1KY,A2KX,A2KY,A3KX,A3KY,A4KX,A4KY,A5KX,A5KY);

    for t = 1:NT-1

        %dsv
        % discountedPayoff = max((S - K * exp(-r * (Tmax - T(t)))), 0);
        % b3x = [discountedPayoff', S'];%In this case rank 2 because of the shape of the condition

        discountedPayoff = max((dsv - K * exp(-r * (Tmax - T(t)))), 0);
        b3x = [discountedPayoff', dsv'];%In this case rank 2 because of the shape of the condition
        
        tol = 1e-5;  % Tolerance for convergence and compression

        [x,y]=CompressData(X,Y,tol);
        % x = X;
        % y = Y;

        % [A,B] = HestonModelOperator(NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);
        
        [AX,AY] = LowRankMatVec(A,B,x,y);
        % [BX,BY] = HestonMatVecBoundaries(NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, K, Tmax, t, T);

        [BX,BY] = HestonMatVecBoundariesLean(b1x,b2x,b3x,b4x,b5x,b1y,b2y,b3y,b4y,b5y);
        %half Euler step
        FX = [(1-r*dt/2)*x,  (dt/2)*AX, dt*BX]; 
        FY = [           y,         AY,    BY];

        %Right hand side vector components
        [BXc,BYc]=CompressData(FX, FY, epsilon);
        %the LHS are operators, not a matrix
        %BXc-->bx
        %BYc-->by
        %?A
        %?B
        % Set GMRES parameters
        % restart = 80;  % Restart after N iterations
        max_iter = iterations;  % Maximum number of iterations

        %%x and y, old values, the initial guesses
                 %GMRES_LowRankV01(x,y, A, B, r, BXc, BYc, x, y, restart, tol, max_iter, dt)
                 %GMRES_LowRankV01(x,y, A, B, r, BX,  BY , x0, y0, restart, tol, max_iter, dt)
        [X, Y] = GMRES_LowRankV01(x, y, A, B, r, BXc, BYc, x, y, restart, tol, max_iter, dt);
        %[X, Y] = GMRES_XYv01(x, y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, K, Tmax, t, T, BXc, BYc, x, y, restart, tol, max_iter, dt);

    end
    [X,Y]=CompressData(X, Y, epsilon);
    U=X*Y';
end


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
        tol = 1e-3;  % Tolerance for convergence and compression

        [x,y]=CompressData(X,Y,tol);
        % x = X;
        % y = Y;

        [AX,AY] = HestonMatVec(x,y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);
        [BX,BY] = HestonMatVecBoundaries(NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, K, Tmax, t, T);

        %half Euler step
        FX = [(1+r*dt/2)*x, (-dt/2)*AX, BX]; 
        FY = [           y,         AY, BY];

        %Right hand side vector components
        [BXc,BYc]=CompressData(FX, FY, epsilon);

        %the LHS are operators, not a matrix


        % Set GMRES parameters
        restart = 3;  % Restart after N iterations (example value)
        
        max_iter = 10;  % Maximum number of iterations
        %%BXc and BYc are the components of the RHS vector
        %%x and y, old values, the initial guesses
        % [x,y]=CompressData(x,y,tol);


               %%GMRES_XYv01(x, y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, K, Tmax, t, T, BX, BY, x0, y0, restart, tol, max_iter) 
                %GMRES_XYv01(x, y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, beta, K, Tmax, t, T, BX, BY, x0, y0, restart, tol, max_iter)               
        [X, Y] = GMRES_XYv01(x, y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, K, Tmax, t, T, BXc, BYc, x, y, restart, tol, max_iter);

    end    
    U=X*Y';
end






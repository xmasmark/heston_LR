function U = HestonExplicitClassicCNXYRC03(params,K,r,q,S,V,T, mode, iterations, restart)

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


        % % % % % % lhs_matrix = ((1+(dt*r/2))*eye(NS*NV) - (dt/2)*AK);  
        % % % % % % rhs_vector = ((1-(dt*r/2))*U_vec + (dt/2)*Av) + dt*Bv;
        % % % % % % 
        % % % % % % % Solve for U_vec at the next time step
        % % % % % % if mode == 0
        % % % % % %     U_vec = (1-dt*r)*U_vec+dt*(Av + Bv);
        % % % % % % end
        % % % % % % if mode == 1
        % % % % % %     U_vec = lhs_matrix \ rhs_vector;
        % % % % % % end
        % % % % % % if mode == 2
        % % % % % %     % Set initial guess for GMRES
        % % % % % %     x0 = U_vec;  % The current solution vector
        % % % % % % 
        % % % % % %     % Set GMRES parameters
        % % % % % %     restart = 80;  % Restart after 20 iterations (example value)
        % % % % % %     tol = 1e-5;  % Tolerance for convergence
        % % % % % %     max_iter = 100;  % Maximum number of iterations
        % % % % % % 
        % % % % % %     % Solve using GMRES
        % % % % % %     warning('off', 'all')
        % % % % % %     U_vec = restarted_gmres(lhs_matrix, rhs_vector, x0, restart, tol, max_iter);
        % % % % % %     warning('on', 'all')
        % % % % % % end



    for t = 1:NT-1
        tol = 1e-5;  % Tolerance for convergence and compression

        [x,y]=CompressData(X,Y,tol);
        % x = X;
        % y = Y;

        [A,B] = HestonModelOperator(NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);
        [AX,AY] = LowRankMatVec(A,B,x,y);

        [AXold,AYold] = HestonMatVec(x,y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);


        [BX,BY] = HestonMatVecBoundaries(NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, K, Tmax, t, T);

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
        [X, Y] = GMRES_LowRankV01(x,y, A, B, r, BXc, BYc, x, y, restart, tol, max_iter, dt);
        %[X, Y] = GMRES_XYv01(x, y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, K, Tmax, t, T, BXc, BYc, x, y, restart, tol, max_iter, dt);

    end    
    U=X*Y';
end

function [X, Y] = ALSOptimization(x, y, FX, FY, BX, BY, epsilon)
% ALS stands for Alternating Linear Scheme
% the concept is to find X and Y solutions as an iterative process
% keeping one dependent variable (V1) fixed at each time and solving the
% linear system on the other variable (V2), then substituting the result
% and iterating on the first (V1)
% the accuracy of the iteration needs to be assessed on the residuals
% I would also add a parameter containing the max number of iterations


end




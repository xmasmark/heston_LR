function U = HestonExplicitClassicCNALSDev01(params,K,r,q,S,V,T, mode, iterations, restart)

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
        tol = 1e-5;  % Tolerance for convergence and compression

        [x,y]=CompressData(X,Y,tol);

        [A,B] = HestonModelOperator(NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);
        [AX,AY] = LowRankMatVec(A,B,x,y);
        
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

        [X, Y] = GMRES_LowRankV01(x,y, A, B, r, BXc, BYc, x, y, restart, tol, max_iter, dt);
        [xALS,yALX]=ALSOptimization(A,B,x,y);
    end    
    U=X*Y';
end

% ALS stands for Alternating Linear Scheme  
% the concept is to find X and Y solutions as an iterative process
% keeping one dependent variable (V1) fixed at each time and solving the
% linear system on the other variable (V2), then substituting the result
% and iterating on the first (V1)
% the accuracy of the iteration needs to be assessed on the residuals
% I would also add a parameter containing the max number of iterations

function [X, Y] = ALSOptimization(A, B, x, y, epsilon)
    
    %solving for X
    YB = LowRankMatVecSingle(B,y);
    YBY=(YB)'*y; %the result is a vector with length R.


end

function [P] = LowRankMatVecSingle(M,LR)

    szA = size(M);

    P = [];

    for n = 1:szA(3)
        P=[P,M(:,:,n)*LR];
    end
end

% function [xl,yl] = LowRankMatVec(A,B,x,y)
% 
%     szA = size(A);
%     szB = size(B);
% 
%     xl = [];
%     yl = [];
% 
%     for n = 1:szA(3)
%         xl=[xl,A(:,:,n)*x];
%     end
% 
%     for n = 1:szB(3)
%         yl=[yl,B(:,:,n)*y];
%     end
% 
%     % 
%     % X1 = A1KX*x; %first slice of A A(:,:,1)
%     % Y1 = A1KY*y; %first slice of B B(:,:,1)
%     % X2 = A2KX*x;
%     % Y2 = A2KY*y;
%     % X3 = A3KX*x;
%     % Y3 = A3KY*y;
%     % X4 = A4KX*x;
%     % Y4 = A4KY*y;
%     % X5 = A5KX*x;
%     % Y5 = A5KY*y;
%     % 
%     % xl = [X1,X2,X3,X4,X5];
%     % yl = [Y1,Y2,Y3,Y4,Y5];
% end



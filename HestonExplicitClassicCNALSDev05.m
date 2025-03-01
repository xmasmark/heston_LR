function U = HestonExplicitClassicCNALSDev05(params,K,r,q,S,V,T, mode, iterations, restart)

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

    d1sM = MDerivativeVM(NS, ds, 0, 1);
    d2sM = MSecondDerivativePlusCVM(NS, ds, 0, 1);
    
    d1vM = MDerivativeVM(NV, dv, 1, 1);
    d2vM = MSecondDerivativePlusCVM(NV, dv, 1, 1);

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

    % SMatrix = diag(S);
    % S2Matrix = diag(S.^2);
    % VMatrix = diag(V);

    b1x = zeros(NS,1);
    b1x(NS)=(1/ds);
    b1y = ones(NV,1);

    b2x = zeros(NS,1);
    b2x(NS)=(1/ds);
    b2y = ones(NV,1);
    
    b3y = zeros(NV,2); %2 columns because this is rank 2
    b3y(1,1)=1;
    b3y(NV,2)=1;
   
    b5x = zeros(NS,1);
    b5x(NS)=(1/ds);
    b5y = ones(NV,1);

    b1x = (r-q) * SMatrix * b1x;
    b2x = 0.5 * S2Matrix * b2x;
    b2y = VMatrix'*b2y;
    b3y = (kappa * (theta * eye(NV) - VMatrix) - lambda * VMatrix)' * b3y;
    b4x = zeros(NS,1);
    b4y = zeros(NV,1);
    
    b4y = 0.5 * sigma^2*VMatrix'*b4y;
    b5x = rho * sigma * SMatrix * b5x;
    b5y = VMatrix'*b5y;

    % xALS = X;
    % yALS = Y;
    x = X;
    y = Y;

    xALS = [X,zeros(NS,6)];
    yALS = [Y,zeros(NV,6)];
    % xALS = X;
    % yALS = Y;

    [A,B] = HestonModelOperator(NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);

    %Adt= each slice of A is multiplied by (dt/2)
    %also I need to add (1-dt*r/2)I(nsxns) and this will be another stack
    %so the rank will be R+1 rather than R
    %Bdt will be a sline I(nvxnv) and the rest stays
    %in both cases, the last slice goes last

    % Ap = APrepared(A,dt,r);
    % Bp = BPrepared(B);
    
    Ap = APrepared(A,dt,r);
    Bp = BPrepared(B);

    for t = 1:NT-1

        [AX,AY] = LowRankMatVec(Ap,Bp,xALS,yALS);
        %BX and BY are constants at each iteration
        [BX,BY] = HestonMatVecBoundaries(NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, K, Tmax, t, T);

        %half Euler step
        FX = [(1-r*dt/2)*xALS,  (dt/2)*AX, dt*BX]; 
        FY = [           yALS,         AY,    BY];

        %Right hand side vector components
        [BXc,BYc]=CompressData(FX, FY, epsilon);

        max_iter = iterations;  % Maximum number of iterations

        [xALS,yALS]=ALSOptimizationW(Ap, Bp, xALS, yALS, BXc, BYc, epsilon, max_iter, restart);



        % % % % % % % tol = 1e-5;  % Tolerance for convergence and compression
        % % % % % % % discountedPayoff = max((S - K * exp(-r * (Tmax - T(t)))), 0);
        % % % % % % % b3x = [discountedPayoff', S'];%In this case rank 2 because of the shape of the condition
        % % % % % % % 
        % % % % % % % %tol = 1e-5;  % Tolerance for convergence and compression
        % % % % % % % [x,y]=CompressData(X,Y,tol);
        % % % % % % % [AX,AY] = LowRankMatVec(A,B,x,y);
        % % % % % % % [BX,BY] = HestonMatVecBoundariesLean(b1x,b2x,b3x,b4x,b5x,b1y,b2y,b3y,b4y,b5y);
        % % % % % % % %half Euler step
        % % % % % % % FX = [(1-r*dt/2)*x,  (dt/2)*AX, dt*BX]; 
        % % % % % % % FY = [           y,         AY,    BY];
        % % % % % % % %Right hand side vector components
        % % % % % % % [BXc,BYc]=CompressData(FX, FY, tol);
        % % % % % % % 
        % % % % % % % max_iter = iterations;  % Maximum number of iterations
        % % % % % % % 
        % % % % % % % [X, Y] = GMRES_LowRankV01(x,y, A, B, r, BXc, BYc, x, y, restart, tol, max_iter, dt);
        % % % % % % % 

    end    
    %U=X*Y';
    U= xALS*yALS';
    % U2 = X*Y';
end


% ALS stands for Alternating Linear Scheme  
% the concept is to find X and Y solutions as an iterative process
% keeping one dependent variable (V1) fixed at each time and solving the
% linear system on the other variable (V2), then substituting the result
% and iterating on the first (V1)
% the accuracy of the iteration needs to be assessed on the residuals
% I would also add a parameter containing the max number of iterations

function [X, Y] = ALSOptimizationW(A, B, x, y, BX, BY, epsilon, max_iter, restart)

    % residual =  ALSEnergyPlus(A, B, x, y, dt, r, BX, BY, epsilon);

    x_opt = x;
    y_opt = y;
    n = 1;

    convergence_iterations = 3;

    residual =  ALSEnergyPlus(A, B, x, y, BX, BY);
    [x_opt, y_opt] = ALSOptimizationV04(A, B, x_opt, y_opt, BX, BY, epsilon, max_iter, restart);

    while abs(residual) > epsilon 
        [x_opt, y_opt] = ALSOptimizationV04(A, B, x_opt, y_opt, BX, BY, epsilon, max_iter, restart);
        residual = ALSEnergyPlus(A, B, x_opt, y_opt, BX, BY);
        n = n+1;

        if n > convergence_iterations
            break
        end
    end

    X = x_opt;
    Y = y_opt;
end

function [Adt] = APrepared(A,dt, r)

    s = size(A);

    Adt(:,:,1)=(-dt/2)*A(:,:,1);
    Adt(:,:,2)=(-dt/2)*A(:,:,2);
    Adt(:,:,3)=(-dt/2)*A(:,:,3);
    Adt(:,:,4)=(-dt/2)*A(:,:,4);
    Adt(:,:,5)=(-dt/2)*A(:,:,5);
    
    Last = (1-dt*r/2)*eye(s(1),s(1));
    Adt(:,:,6)=Last;
end

function [Bdt] = BPrepared(B)

    s = size(B);

    Bdt(:,:,1)=B(:,:,1);
    Bdt(:,:,2)=B(:,:,2);
    Bdt(:,:,3)=B(:,:,3);
    Bdt(:,:,4)=B(:,:,4);
    Bdt(:,:,5)=B(:,:,5);
    
    Last = eye(s(1),s(1));
    Bdt(:,:,6)=Last;
end

function residual = ALSEnergyPlus(A, B, x, y, BXc, BYc)

    szA = size(A);
    R = szA(3);
    sx=size(x);
    r = sx(2);

    YB = pagemtimes(B,y);
    YBt = permute(YB,[2,1,3]); %transpose each layer of YB
    YBY = pagemtimes(YBt,y);

    XA = pagemtimes(A,x);
    XAt = permute(XA,[2,1,3]); %transpose each layer of XA
    XAX = pagemtimes(XAt,x);

    xv = reshape(YBY,r*r*R,1);
    yv = reshape(XAX,r*r*R,1);

    one_side = xv'*yv;
    
    xBx = x'*BXc;
    yBYc = y'*BYc;

    s= size(xBx);
    one = reshape(xBx,s(1)*s(2),1);
    two = reshape(yBYc,1,s(1)*s(2));
    other_side = two*one;

    residual = 0.5*(one_side) - other_side;
end

function [X, Y] = ALSOptimizationV04(A, B, x, y, BXc, BYc, epsilon, max_iter, restart)

    szA = size(A);
    NS = szA(1);
    R = szA(3);
    sx=size(x);
    r = sx(2);

    sizeB=size(B);
    NV = sizeB(1);

    %solving for X
    %left hand side part of CN also called Y*B*Y in the notes
    YB = pagemtimes(B,y);
    YBt = permute(YB,[2,1,3]);
    YBY = pagemtimes(YBt,y);

    AR = reshape(A,NS*NS,R);
    YBYR = reshape(YBY,r*r,R);

    ah = AR*YBYR';
    A_hat = reshape(ah,NS,NS,r,r);

    Y_BYc = y'*BYc;
    b_hat = BXc*Y_BYc';

    %given the values I get, the following step is definitely wrong
    A_hat_matrix = A_hat;
    b_hat_vector = b_hat;

    if r >= 2
        %A_hat -- sizes are NS, NS, r and r
        A_hatP = permute(A_hat,[1,3,2,4]);
        A_hat_matrix = reshape(A_hatP,NS*r,NS*r);
        b_hat_vector = reshape(b_hat,NS*r,1);
    end

    x0=reshape(x,NS*r,1);

    %[x, flag, relres, iter]
    %[X_Opt, flag, relres, iter] = gmres_simple(A_hat_matrix, b_hat_vector, epsilon, max_iter);
    X_Opt = restarted_gmres(A_hat_matrix, b_hat_vector, x0, restart, epsilon, max_iter);
    % X_Opt = A_hat_matrix \ b_hat_vector;

    X_OptR = reshape(X_Opt,NS,r);


    % oldEnergy = ALSEnergyPlus(A, B, x, y, BXc, BYc);
    % newEnergy = ALSEnergyPlus(A, B, X_OptR, y, BXc, BYc);

    %solving for Y
    XA = pagemtimes(A,X_OptR);
    XAt = permute(XA,[2,1,3]);
    XAX = pagemtimes(XAt, X_OptR);

    BR = reshape(B,NV*NV,R);
    XAXR = reshape(XAX,r*r,R);

    ahY = BR*XAXR';
    A_hatY = reshape(ahY,NV,NV,r,r);

    X_BXc=X_OptR'*BXc;
    b_hatY_vector = X_BXc*BYc';

    if r>1
        %bHat=X_OptR'*BXc*BYc';
        %A_hatY -- sizes are NV, NV, r and r
        A_hatYP = permute(A_hatY,[1,3,2,4]);
        %A_hatYP = permute(A_hatY,[2,4,1,3]);
        A_hatY_matrix = reshape(A_hatYP,NV*r,NV*r);
        b_hatY_vector = reshape(b_hatY_vector,NV*r,1);
        %Y_Opt = A_hatY_matrix \ b_hatY_vector;
    end
    if r==1
       %A_hatY_matrix = reshape(A_hatY,NV*r,NV*r);
       A_hatY_matrix = A_hatY;
       b_hatY_vector = b_hatY_vector';
       %Y_Opt = A_hatY \ b_hatY_vector';
    end

    y0 = reshape(y,NV*r,1);

    % Y_Opt = A_hatY_matrix \ b_hatY_vector;
    %[X_Opt, flag, relres, iter] 
    %[Y_Opt, flag, relres, iter] = gmres_simple(A_hatY_matrix, b_hatY_vector, epsilon, max_iter);
    Y_Opt = restarted_gmres(A_hatY_matrix, b_hatY_vector, y0, restart, epsilon, max_iter);
   
    X_Opt = reshape(X_Opt,NS,r);
    Y_Opt = reshape(Y_Opt,NV,r);

    %calculation of residuals to control the progress

    X=X_Opt;
    Y=Y_Opt;

end



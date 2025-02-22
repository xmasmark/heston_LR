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
    
    % SMatrix = diag(S);
    % S2Matrix = diag(S.^2);
    % VMatrix = diag(V);

    xALS = X;
    yALS = Y;

    for t = 1:NT-1
        tol = 1e-5;  % Tolerance for convergence and compression

        [x,y]=CompressData(X,Y,tol);
        [xALS,yALS]=CompressData(xALS,yALS,tol);

        [A,B] = HestonModelOperator(NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);
        [AX,AY] = LowRankMatVec(A,B,x,y);
        [AXALS,AYALS] = LowRankMatVec(A,B,xALS,yALS);
        
        [BX,BY] = HestonMatVecBoundaries(NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, K, Tmax, t, T);

        %half Euler step
        FX = [(1-r*dt/2)*x,  (dt/2)*AX, dt*BX]; 
        FY = [           y,         AY,    BY];

        FXALS = [(1-r*dt/2)*xALS,  (dt/2)*AXALS, dt*BX]; 
        FYALS = [           yALS,         AYALS,    BY];
        
        %Right hand side vector components
        [BXc,BYc]=CompressData(FX, FY, epsilon);
        [BXcALS,BYcALS]=CompressData(FXALS, FYALS, epsilon);

        %the LHS are operators, not a matrix
        %BXc-->bx
        %BYc-->by
        %?A
        %?B
        % Set GMRES parameters
        % restart = 80;  % Restart after N iterations
        max_iter = iterations;  % Maximum number of iterations

        [X, Y] = GMRES_LowRankV01(x,y, A, B, r, BXc, BYc, x, y, restart, tol, max_iter, dt);

        [xALS,yALS]=ALSOptimizationW(A, B, x, y, BXc, BYc, epsilon,max_iter);

        
    end    
    %U=X*Y';
    U= xALS*yALS';
end

% ALS stands for Alternating Linear Scheme  
% the concept is to find X and Y solutions as an iterative process
% keeping one dependent variable (V1) fixed at each time and solving the
% linear system on the other variable (V2), then substituting the result
% and iterating on the first (V1)
% the accuracy of the iteration needs to be assessed on the residuals
% I would also add a parameter containing the max number of iterations

function [X, Y] = ALSOptimizationW(A, B, x, y, BXc, BYc, epsilon, max_iter)

    residual = ALSResidualCalculation(A, B, x, y, BXc, BYc, epsilon);

    x_opt = x;
    y_opt = y;
    n = 1;

    restart = 5;

    %ALSOptimizationV04(A, B, x, y, BXc, BYc, epsilon, max_iter, restart)
    [x_opt, y_opt] = ALSOptimizationV04(A, B, x_opt, y_opt, BXc, BYc, epsilon, max_iter, restart);

    while abs(residual) > epsilon 
        [x_opt, y_opt] = ALSOptimizationV04(A, B, x_opt, y_opt, BXc, BYc, epsilon, max_iter, restart);
        residual = ALSResidualCalculation(A, B, x_opt, y_opt, BXc, BYc, epsilon);
        n = n+1;
        if n>max_iter
            break
        end
    end

    X=x_opt;
    Y = y_opt;
    %[X, Y] = ALSOptimizationV02(A, B, x, y, BXc, BYc, epsilon);
end

function residual = ALSResidualCalculation(A, B, x, y, BXc, BYc, epsilon)

    szA = size(A);
    R = szA(3);
    sx=size(x);
    r = sx(2);

    YB =LowRankMatVecStacked(B,y);
    YBY =LowRankMatVecStacked(YB,y);
    XA =LowRankMatVecStacked(A,x);
    XAX =LowRankMatVecStacked(XA,x);

    xv = reshape(YBY,r*r*R,1);
    yv = reshape(XAX,r*r*R,1);

    one_side = xv'*yv;
    
    xBx = x'*BXc;
    yBYc = y'*BYc;

    s= size(xBx);
    one = reshape(xBx,s(1)*s(2),1);
    two = reshape(yBYc,1,s(1)*s(2));
    other_side = two*one;

    residual = one_side - other_side;
end


function [X, Y] = ALSOptimizationV04(A, B, x, y, BXc, BYc, epsilon, max_iter, restart)

    szA = size(A);
    NS = szA(1);
    R = szA(3);
    sx=size(x);
    r = sx(2);

    sizeB=size(B);
    NV = sizeB(1);

    %sol    ving for X
    %left hand side part of CN also called Y*B*Y in the notes
    YB =LowRankMatVecStacked(B,y);
    %ybyp = permute(y' * reshape(B, NV, []), [1, 3, 2]);
    %YBpp = reshape(reshape(B, NV, NV*R) * y', NV, r, R);
    %P(:,:,n)=(LR'*M(:,:,n))';

    YBY =LowRankMatVecStacked(YB,y);

    %AR = reshape(A,szA(1)*szA(1),szA(3));
    AR = reshape(A,NS*NS,R);
    %szB =size(YBY);
    YBYR = reshape(YBY,r*r,R);

    ah = AR*YBYR';
    A_hat = reshape(ah,NS,NS,r,r);

    Y_BYc = y'*BYc;
    b_hat = BXc*Y_BYc';

    %given the values I get, the following step is definitely wrong
    A_hat_matrix = A_hat;
    b_hat_vector = b_hat;

    % if r == 1
    %     A_hat_matrix = A_hat;
    %     b_hat_vector = b_hat;
    % end
    if r > 2
        %A_hat -- sizes are NS, NS, r and r
        A_hatP = permute(A_hat,[1,3,2,4]);
        A_hat_matrix = reshape(A_hatP,NS*r,NS*r);
        b_hat_vector = reshape(b_hat,NS*r,1);
    end

    %[x, flag, relres, iter]
    %[X_Opt, flag, relres, iter] = gmres_simple(A_hat_matrix, b_hat_vector, epsilon, max_iter);
    X_Opt = restarted_gmres(A_hat_matrix, b_hat_vector, x, restart, epsilon, max_iter);
    %X_Opt = A_hat_matrix \ b_hat_vector;

    %solving for Y
    if r == 1
        X_OptR = X_Opt;
    end
    if r > 1
        X_OptR = reshape(X_Opt,NS,r);
    end

    XA =LowRankMatVecStacked(A,X_OptR);
    XAX =LowRankMatVecStacked(XA,X_OptR);
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

    %Y_Opt = A_hatY_matrix \ b_hatY_vector;
    %[X_Opt, flag, relres, iter] 
    %[Y_Opt, flag, relres, iter] = gmres_simple(A_hatY_matrix, b_hatY_vector, epsilon, max_iter);
    Y_Opt = restarted_gmres(A_hatY_matrix, b_hatY_vector, y, restart, epsilon, max_iter);
   
    X_Opt = reshape(X_Opt,NS,r);
    Y_Opt = reshape(Y_Opt,NV,r);

    %calculation of residuals to control the progress

    X=X_Opt;
    Y=Y_Opt;

end


function [X, Y] = ALSOptimizationV03(A, B, x, y, BXc, BYc, epsilon, max_iter)
    

    szA = size(A);
    NS = szA(1);
    R = szA(3);
    sx=size(x);
    r = sx(2);

    sizeB=size(B);
    NV = sizeB(1);

    %sol    ving for X
    %left hand side part of CN also called Y*B*Y in the notes
    YB =LowRankMatVecStacked(B,y);
    %ybyp = permute(y' * reshape(B, NV, []), [1, 3, 2]);
    %YBpp = reshape(reshape(B, NV, NV*R) * y', NV, r, R);
    %P(:,:,n)=(LR'*M(:,:,n))';

    YBY =LowRankMatVecStacked(YB,y);

    %AR = reshape(A,szA(1)*szA(1),szA(3));
    AR = reshape(A,NS*NS,R);
    %szB =size(YBY);
    YBYR = reshape(YBY,r*r,R);

    ah = AR*YBYR';
    A_hat = reshape(ah,NS,NS,r,r);

    Y_BYc = y'*BYc;
    b_hat = BXc*Y_BYc';

    %given the values I get, the following step is definitely wrong
    A_hat_matrix = A_hat;
    b_hat_vector = b_hat;

    % if r == 1
    %     A_hat_matrix = A_hat;
    %     b_hat_vector = b_hat;
    % end
    if r > 2
        %A_hat -- sizes are NS, NS, r and r
        A_hatP = permute(A_hat,[1,3,2,4]);
        A_hat_matrix = reshape(A_hatP,NS*r,NS*r);
        b_hat_vector = reshape(b_hat,NS*r,1);
    end

    %[x, flag, relres, iter]
    [X_Opt, flag, relres, iter] = gmres_simple(A_hat_matrix, b_hat_vector, epsilon, max_iter);
    %X_Opt = A_hat_matrix \ b_hat_vector;

    %solving for Y
    if r == 1
        X_OptR = X_Opt;
    end
    if r > 1
        X_OptR = reshape(X_Opt,NS,r);
    end

    XA =LowRankMatVecStacked(A,X_OptR);
    XAX =LowRankMatVecStacked(XA,X_OptR);
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

    %Y_Opt = A_hatY_matrix \ b_hatY_vector;
    %[X_Opt, flag, relres, iter] 
    [Y_Opt, flag, relres, iter] = gmres_simple(A_hatY_matrix, b_hatY_vector, epsilon, max_iter);
   
    X_Opt = reshape(X_Opt,NS,r);
    Y_Opt = reshape(Y_Opt,NV,r);

    %calculation of residuals to control the progress

    X=X_Opt;
    Y=Y_Opt;

end


function [X, Y] = ALSOptimizationV02(A, B, x, y, BXc, BYc, epsilon, max_iter)
    
    %solving for X
    %left hand side part of CN also called Y*B*Y in the notes
    YB =LowRankMatVecStacked(B,y);
    YBY =LowRankMatVecStacked(YB,y);
    %right hand side of CN Y*BY*BX
    %YBY=(YB)'*y; %the result is a vector with length R.

    % % % OK, A(Ns,Ns,R)  YBY(r,r,R)
    % % % A = reshape(A,[Ns*Ns, R]);
    % % % YBY = reshape(YBY,[r*r,R]);
    % % % Ahat = A*YBY';
    % % % Ahat = reshape(Ahat,[Ns,Ns,r,r]);

    %change the sizes into NS, NV, r and R
    %remove subroutines

    szA = size(A);
    NS = szA(1);
    R = szA(3);
    sx=size(x);
    r = sx(2);

    sizeB=size(B);
    NV = sizeB(1);

    %AR = reshape(A,szA(1)*szA(1),szA(3));
    AR = reshape(A,NS*NS,R);
    %szB =size(YBY);
    YBYR = reshape(YBY,r*r,R);

    ah = AR*YBYR';
    A_hat = reshape(ah,NS,NS,r,r);

    Y_BYc = y'*BYc;
    b_hat = BXc*Y_BYc';

    %given the values I get, the following step is definitely wrong
    A_hat_matrix = A_hat;
    b_hat_vector = b_hat;

    % if r == 1
    %     A_hat_matrix = A_hat;
    %     b_hat_vector = b_hat;
    % end
    if r > 2
        %A_hat -- sizes are NS, NS, r and r
        A_hatP = permute(A_hat,[1,3,2,4]);
        A_hat_matrix = reshape(A_hatP,NS*r,NS*r);
        b_hat_vector = reshape(b_hat,NS*r,1);
    end

    %[x, flag, relres, iter]
    [X_Opt, flag, relres, iter] = gmres_simple(A_hat_matrix, b_hat_vector, epsilon, max_iter);
    %X_Opt = A_hat_matrix \ b_hat_vector;

    %solving for Y
    if r == 1
        X_OptR = X_Opt;
    end
    if r > 1
        X_OptR = reshape(X_Opt,NS,r);
    end

    XA =LowRankMatVecStacked(A,X_OptR);
    XAX =LowRankMatVecStacked(XA,X_OptR);
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

    %Y_Opt = A_hatY_matrix \ b_hatY_vector;
    %[X_Opt, flag, relres, iter] 
    [Y_Opt, flag, relres, iter] = gmres_simple(A_hatY_matrix, b_hatY_vector, epsilon, max_iter);
   
    X_Opt = reshape(X_Opt,NS,r);
    Y_Opt = reshape(Y_Opt,NV,r);

    %calculation of residuals to control the progress

    X=X_Opt;
    Y=Y_Opt;

end




function [X, Y] = ALSOptimization(A, B, x, y, BXc, BYc, epsilon)
    
    %solving for X
    %left hand side part of CN also called Y*B*Y in the notes
    YB =LowRankMatVecStacked(B,y);
    YBY =LowRankMatVecStacked(YB,y);
    %right hand side of CN Y*BY*BX
    %YBY=(YB)'*y; %the result is a vector with length R.

    % % % OK, A(Ns,Ns,R)  YBY(r,r,R)
    % % % A = reshape(A,[Ns*Ns, R]);
    % % % YBY = reshape(YBY,[r*r,R]);
    % % % Ahat = A*YBY';
    % % % Ahat = reshape(Ahat,[Ns,Ns,r,r]);

    %change the sizes into NS, NV, r and R
    %remove subroutines

    szA = size(A);
    NS = szA(1);
    R = szA(3);
    sx=size(x);
    r = sx(2);

    sizeB=size(B);
    NV = sizeB(1);

    %AR = reshape(A,szA(1)*szA(1),szA(3));
    AR = reshape(A,NS*NS,R);
    szB =size(YBY);
    %YBYR = reshape(YBY,szB(1)*szB(2),szB(3));
    YBYR = reshape(YBY,r*r,R);

    ah = AR*YBYR';
    %szah = size(ah);
    %A_hat = reshape(ah,szA(1),szA(1),r,r);
    A_hat = reshape(ah,NS,NS,r,r);

    Y_BYc = y'*BYc;
    b_hat = BXc*Y_BYc';

    ahs = size(A_hat);

    dim = ndims(A_hat);

    %given the values I get, the following step is definitely wrong
    if dim < 3
        A_hat_matrix = A_hat;
        b_hat_vector = b_hat;
    end
    
    if dim > 3
        %A_hat_matrix = reshape(A_hat,ahs(1)*ahs(3),ahs(1)*ahs(3));
        A_hatP = permute(A_hat,[1,3,2,4]);
        %A_hat_matrix = reshape(A_hatP,ahs(1)*ahs(3),ahs(1)*ahs(3));
        A_hat_matrix = reshape(A_hatP,NS*r,NS*r);
        %b_hat_vector = reshape(b_hat,ahs(1)*ahs(3),1);
        b_hat_vector = reshape(b_hat,NS*r,1);
    end
   
    X_Opt = A_hat_matrix \ b_hat_vector;

    %solving for Y
    if dim < 3
        X_OptR = X_Opt;
    end
    if dim > 3
        X_OptR = reshape(X_Opt,ahs(1),ahs(3));
    end

    % AR = reshape(A,szA(1)*szA(1),szA(3));
    % szB =size(YBY);
    % YBYR = reshape(YBY,szB(1)*szB(2),szB(3));


    XA =LowRankMatVecStacked(A,X_OptR);
    XAX =LowRankMatVecStacked(XA,X_OptR);
    BR = reshape(B,NV*NV,R);
    %YB =LowRankMatVecStacked(B,y);
    %PROBLEM HERE: I need to connect XAX and BY and then the AHat matrix is
    %calculated
    %ah = AR*YBYR';
    XAXR = reshape(XAX,r*r,R);
    %YBR = reshape(YB,NV*r,R);

    %ah = AR*YBYR';
    ahY = BR*XAXR';
    A_hatY = reshape(ahY,NV,NV,r,r);

    % Y_BYc = y'*BYc;
    % b_hat = BXc*Y_BYc';
    X_BXc=X_OptR'*BXc;
    b_hatY_vector = X_BXc*BYc';

    if r>1
        bHat=X_OptR'*BXc*BYc';
        A_hatYP = permute(A_hatY,[1,3,2,4]);
        A_hatY_matrix = reshape(A_hatYP,NV*r,NV*r);
        b_hatY_vector = reshape(b_hatY_vector,NV*r,1);
        Y_Opt = A_hatY_matrix \ b_hatY_vector;
    end
    if r==1
       A_hatY_matrix = reshape(A_hatY,NV*r,NV*r);
       Y_Opt = A_hatY \ b_hatY_vector';
    end


    % % % if r > 1
    % % %     Y_Opt = reshape(Y_Opt,NV,r);
    % % % end
    
    X_Opt = reshape(X_Opt,NS,r);
    Y_Opt = reshape(Y_Opt,NV,r);

    % % % % % % % 
    % % % % % % % %reshape
    % % % % % % % szXAX =size(XAX);
    % % % % % % % XAXR = reshape(XAX,szXAX(1)*szXAX(2),szXAX(3));
    % % % % % % % szY_optB=size(YB);
    % % % % % % % Y_optBR=reshape(YB,szY_optB(1)*szY_optB(2),szY_optB(3));
    % % % % % % % 
    % % % % % % % % XA =LowRankMatVecStacked(A,x);
    % % % % % % % % XAX =LowRankMatVecStacked(XA,x);
    % % % % % % % % Y_optB =LowRankMatVecStacked(B,Y_Opt);
    % % % % % % % % %reshape
    % % % % % % % % szXAX =size(XAX);
    % % % % % % % % XAXR = reshape(XAX,szXAX(1)*szXAX(2),szXAX(3));
    % % % % % % % % szY_optB=size(Y_optB);
    % % % % % % % % Y_optBR=reshape(Y_optB,szY_optB(1)*szY_optB(2),szY_optB(3));
    % % % % % % % 
    % % % % % % % %create the new A_hat and b_hat...
    % % % % % % % ahx=XAXR*Y_optBR';
    % % % % % % % 
    % % % % % % % bhatx=x'*BXc*BYc';

    %calculation of residuals to control the progress

    X=X_Opt;
    Y=Y_Opt;

end


%this should return a 3d structure
%on one side there is the operator M, on the other side the low rank
%structure LR
%the operator M is composed by R layers, one layer for each
%operator component so the result of this should be R new layers where the
%operator has worked each of its components on the low rank structure
function [P] = LowRankMatVecStacked(M,LR)
    szA = size(M);
    P = [];
    for n = 1:szA(3)
        P(:,:,n)=(LR'*M(:,:,n))';
    end
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



function MU = HestonExplicitLowRankRTildaPlusGMRES(params,K,r,q,S,V,T,epsilon)




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
    
    %testMatrix = U-X*Y';
    
    %caching constant elements
    xzeros = zeros(NS,1);
    xzeros2 = zeros(NS,1);
    yzeros = zeros(NV,1);
    yzeros2 = zeros(NV,1);
    yzeros(NV)=1;
    yzeros2(1)=1;

    for t=1:NT-1

        %the following will be now commented out. The reason for this is
        %the fact that I need to apply the text book bounday condition, so
        %either Dirichlet or Neumann.
        
        %I leave the following two statements
        x = [X];
        y = [Y];

        %%%%%Boundary Conditions start

        % % % % % % % x = [X xzeros2 xzeros];
        % % % % % % % y = [Y yzeros yzeros2];
        % % % % % % % 
        % % % % % % % %first step is to 
        % % % % % % % 
        % % % % % % % sz = size(x);
        % % % % % % % records = sz(1);
        % % % % % % % fields = sz(2);
        % % % % % % % 
        % % % % % % % %U(s,NV) = S(s);
        % % % % % % % 
        % % % % % % % %U(s,1) = ...;
        % % % % % % % %boundary condition for the first column of the matrix U
        % % % % % % % %once this part of the code will reach maturity
        % % % % % % % %both this loop and the previous one will be merged into one
        % % % % % % % 
        % % % % % % % yc=Y(1,:);
        % % % % % % % yc2=Y(2,:);
        % % % % % % % ynv = Y(NV,:);
        % % % % % % % 
        % % % % % % % for j = 2:(records-1)
        % % % % % % %     xc=X(j,:);
        % % % % % % %     xcp1=X(j+1,:);
        % % % % % % %     xcm1=X(j-1,:);
        % % % % % % % 
        % % % % % % %     us1 = xc*yc';
        % % % % % % %     usp1= xcp1*yc';
        % % % % % % %     us2 = xc*yc2';
        % % % % % % %     usm1= xcm1*yc';
        % % % % % % % 
        % % % % % % %     %new value, hence the n prefix:
        % % % % % % %     nus1=us1*(1-r*dt-kappa*theta*dt/dv)+dt*(0.5*(r-q)*j*usp1-usm1)+kappa*theta*dt/dv*us2;
        % % % % % % % 
        % % % % % % %     %in acc I need to store the cross 
        % % % % % % %     %product between X and the first column of Y'
        % % % % % % % 
        % % % % % % %     acc = xc*yc';
        % % % % % % % 
        % % % % % % %     x(j,fields)=nus1-acc;
        % % % % % % % 
        % % % % % % %     acc2 = xc*ynv';            
        % % % % % % %     x(j,fields-1)=S(j)-acc2;
        % % % % % % % end
        % % % % % % % 
        % % % % % % % %Last row of X that couldn't be covered in the previous loop
        % % % % % % % xc=X(records,:);
        % % % % % % % acc2 = xc*ynv';            
        % % % % % % % x(records,fields-1)=S(records)-acc2;
        % % % % % % % 
        % % % % % % % % for i = 1:fields-1
        % % % % % % % %     x(1,i)=0;
        % % % % % % % % end
        % % % % % % % 
        % % % % % % % x(1,:)=0;

        %%%%%Boundary Conditions end

        sz = size(y);
        y_rows = sz(1);
        y_cols = sz(2);

        szx = size(x);
        x_rows = szx(1);
        x_cols = szx(2);

        %Dirichlet boundary condition - correction factor
        %u_left = (u2-u0); -- in matlab the counting starts from 1 so there
        %is no "zero" index, so what in the write up is 0-2 becomes 1-3
        %I don't have u anymore but X and Y, these are two matrices and
        %they require their very own correction factors

        %Boundary Type: 0 --> Dirichlet, 1 --> Neumann
        boundaryType = 0;

        x3 = x(3,:);
        x1 = x(1,:);
        x2 = x(2,:);

        x_last = x(x_rows,:);
        x_lastm2 = x(x_rows-2,:);
        x_lastm1 = x(x_rows-1,:);

        y3 = y(3,:);
        y1 = y(1,:);
        y2 = y(2,:);

        y_last = y(y_rows,:);
        y_lastm2 = y(y_rows-2,:);
        y_lastm1 = y(y_rows-1,:);

        xLeft = 0;
        xRight = 0;
        yLeft = 0;
        yRight = 0;

        if(boundaryType==0)
            %Dirichelet boundary conditions
            xLeft = -(x3-x1);          %this is the left value
            xRight = x_last-x_lastm2;  %this is the right value
    
            yLeft = -(y3-y1);          %this is the left value
            yRight = y_last-y_lastm2;  %this is the right value

        elseif(boundaryType==1)
            %Neumann boundary conditions
            xLeft = (x2-x1)/ds;
            xRight = (x_last-x_lastm1)/ds;

            yLeft = (y2-y1)/dv;
            yRight = (y_last-y_lastm1)/dv;
        end

        %now create the matrix
        xcorrection = zeros(x_rows,x_cols);
        xcorrection(1,:)=xLeft;
        xcorrection(x_rows,:)=xRight;

        %now create the matrix
        ycorrection = zeros(y_rows,y_cols);
        ycorrection(1,:)=yLeft;
        ycorrection(y_rows,:)=yRight;

        %MDerivative(matrix, delta, boundaryConditionLeft, boundaryConditionRight, boundaryLeft, boundaryRight)
        fdX = MDerivative(x,ds,boundaryType,boundaryType,xcorrection);     %first derivative of X
        fdY = MDerivative(y,dv,boundaryType,boundaryType,ycorrection);     %first derivative of Y
        sdX = MSecondDerivative(x,ds);  %second derivative of X
        sdY = MSecondDerivative(y,dv);  %second derivative of Y
    
        spVectorSquare = S'.^2;

        %PWProduct needs to take into account that I am
        %now dealing with X and Y as matrixes and not vectors
        
        %s2sdX=PWProduct(spVectorSquare,sdX);

        s2sdX = spVectorSquare.*sdX;

        vY=V'.*y;

        vsdY=V'.*sdY;

        vfdY=V'.*fdY;

        sfdX=S'.*fdX;

        thetaKappaElement = (kappa*(theta-V')).*fdY;
        
        %firstElement = 0.5*s2sdX*vY';
        X1=0.5*s2sdX;
        Y1=vY;

        %secondElement = 0.5*sigma^2*x*vsdY';
        X2=0.5*sigma^2*x;
        Y2=vsdY;

        % %thirdElement = sigma*sfdX*vfdY';
        X3 = sigma*sfdX; %*0.25; %added 0.25 because excluded from direct formula
        Y3 = vfdY;
        
        %fourthElement = (r-q)*sfdX*Y';
        X4 = (r-q)*sfdX;
        Y4 = y;

        X5 = x;
        Y5 = thetaKappaElement;

        X6 = -r*x;
        Y6 = y;

        dtsrt = sqrt(dt);

        XC = dtsrt*[X1,X2,X3,X4,X5,X6];
        YC = dtsrt*[Y1,Y2,Y3,Y4,Y5,Y6];

        %where should I plug in the call to restarted_gmres_low_rank?
        %what is A? what is B?
    	%u = U;
    	%b = B*u;
        %x0 = zeros(size(b));
        %uM = A\B*u;
        %[U] = restarted_gmres(A, b, x0, 15, epsilon, iterations);

        newX = [x XC];
        newY = [y YC];

        [X,Y]=CompressData(newX,newY,epsilon);

        % try
        %     [X,Y]=CompressData(newX,newY,epsilon);
        % catch
        %     a = 1;
        % end
        
    end

    MU=X*Y';
end

function [X, Y]=CompressData(XC,YC, epsilon)

    %[Qx,Rx]=qr(XC,0,"econ");
    %[Qy,Ry]=qr(YC,0,"econ");

    [Qx,Rx]=qr(XC,"econ");
    [Qy,Ry]=qr(YC,"econ");

    Ri = Rx*Ry';
    
    [rank, u, s, v]=SVDRankPlus(Ri,epsilon);

    %the prefix d stands for "dominant"
    du = u(:,1:rank);
    ds = s(1:rank,1:rank);
    %dv = v(1:rank,:);   
    dv = v(:,1:rank);   

    %creation of compressed matrices
    X=Qx*du*sqrt(ds);
    Y=Qy*dv*sqrt(ds);

end



function D = MDerivative_old(L, delta)
    % Calculate the size of the input matrix L
    [records, fields] = size(L);

    % Initialize the output matrix D with zeros
    D = zeros(records, fields);

    % Calculate the first derivative using diff and vectorized operations
    if records > 2
        D(2:end-1, :) = (L(3:end, :) - L(1:end-2, :)) / (2 * delta);
    end
end


function D = MDerivative_with_Dmitry(matrix, delta, boundaryConditionLeft, boundaryConditionRight, boundaryLeft, boundaryRight)
    %boundaryCondition -- Dirichlet = 0, Neumann = 1
    %boundaryLeft and boundaryRight are the actual data 
    %to insert in the correction vector

    % Calculate the size of the input matrix L
    [records, fields] = size(matrix);
   
    n = records;
    h = delta;

    boundaryVector = zeros(n,1);

    % Create the "Laplacian" matrix
    temporary = zeros(n, n);
    temporary = temporary - diag(ones(n-1, 1), -1);
    temporary = temporary + diag(ones(n-1, 1), 1);
    
    % % Handle boundary conditions
    % temporary(1, 1) = -1;
    % temporary(1, 2) = 1;
    % temporary(n, n-1) = -1;
    % temporary(n, n) = 1;

    if(boundaryConditionLeft == 0)
        temporary(1,1)=0;
        %boundaryVector(1)=-boundaryLeft/(2*h);
        boundaryVector(1)=-boundaryLeft;
    elseif(boundaryConditionLeft == 1)
        temporary(1,1)=-1;
        %boundaryVector(1)=boundaryLeft/(2);
        boundaryVector(1)=boundaryLeft;
    end

    if(boundaryConditionRight == 0)
        temporary(n,n)=0;
        %boundaryVector(n)=boundaryRight/(2*h);
        boundaryVector(n)=boundaryRight;
    elseif(boundaryConditionLeft == 1)
        temporary(n,n)=1;
        %boundaryVector(n)=boundaryRight/(2);
        boundaryVector(n)=boundaryRight;
    end

    % Scale by 1/(2*h) to account for the step size
    temporary = temporary / (2 * h);

    % Compute the derivative
    D = temporary * matrix+boundaryVector;

end

%for now it is better to pre-compute the boundaries in the code above and
%to drive this method through the precomputed boundaries
function D = MDerivative(matrix, delta, boundaryConditionLeft, boundaryConditionRight, boundaryMatrix)
    %boundaryCondition -- Dirichlet = 0, Neumann = 1
    %boundaryLeft and boundaryRight are the actual data 
    %to insert in the correction vector

    % Calculate the size of the input matrix L
    [records, fields] = size(matrix);
   
    n = records;
    h = delta;

    %boundaryVector = zeros(n,1);

    % Create the "Laplacian" matrix
    temporary = zeros(n, n);
    temporary = temporary - diag(ones(n-1, 1), -1);
    temporary = temporary + diag(ones(n-1, 1), 1);
    
    % % Handle boundary conditions
    % temporary(1, 1) = -1;
    % temporary(1, 2) = 1;
    % temporary(n, n-1) = -1;
    % temporary(n, n) = 1;

    if(boundaryConditionLeft == 0)
        temporary(1,1)=0;
        %boundaryVector(1)=-boundaryLeft/(2*h);
        %boundaryVector(1)=-boundaryLeft;
    elseif(boundaryConditionLeft == 1)
        temporary(1,1)=-1;
        %boundaryVector(1)=boundaryLeft/(2);
        %boundaryVector(1)=boundaryLeft;
    end

    if(boundaryConditionRight == 0)
        temporary(n,n)=0;
        %boundaryVector(n)=boundaryRight/(2*h);
        %boundaryVector(n)=boundaryRight;
    elseif(boundaryConditionLeft == 1)
        temporary(n,n)=1;
        %boundaryVector(n)=boundaryRight/(2);
        %boundaryVector(n)=boundaryRight;
    end

    % Scale by 1/(2*h) to account for the step size
    temporary = temporary / (2 * h);

    % Compute the derivative
    D = (temporary * matrix) + boundaryMatrix;

end

%derS=(u(s+1,v)-u(s-1,v))/(2*ds);
%derS2=(u(s+1,v)-2*u(s,v)+u(s-1,v))/(ds^2);

function D = MSecondDerivative_old(L, delta)
    % Calculate the size of the input matrix L
    [records, fields] = size(L);

    % Initialize the output matrix D with zeros
    D = zeros(records, fields);

    % Calculate the second derivative using diff and vectorized operations
    if records > 2
        D(2:end-1, :) = (L(3:end, :) + L(1:end-2, :) - 2 * L(2:end-1, :)) / (delta^2);
    end
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

function [R, uc, sc, vs] = SVDRankPlus(U, epsilon)

    %builtin_rank = rank(U,epsilon);

    [uc,sc,vs]=svd(U, "econ");

    sigSum = 0;

    max=size(sc,2);
    if(size(sc,1)<size(sc,2))
        max = size(sc,1);
    end

    for i = 1: max
        sigSum = sigSum + power(sc(i,i),2);
    end

    threshold = power(epsilon,2)*sigSum;

    rank_home_made = 0;

    errorSquare = 0;

    for i=size(sc,2):-1:1
        errorSquare = errorSquare + power(sc(i,i),2);
        if errorSquare > threshold
            rank_home_made = i + 2;
            break;
        end
    end

    if(rank_home_made>max)
        rank_home_made = max;
    end

    R=rank_home_made;

end

function [xl,yl] = HestonMatVec(x,y)

    %make sure this works -- pass additional parameters like ds, dt, dv,
    %NS, NV, S, V, K, Tmax, T...
    %also, create independent functions for each method call...
    %stand alone functions in separate files

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

    xl = [X1,X2,X3,X4,X5];
    yl = [Y1,Y2,Y3,Y4,Y5];
end
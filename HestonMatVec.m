function [xl,yl] = HestonMatVec(x,y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho)
    %make sure this works -- pass additional parameters like ds, dt, dv,
    %NS, NV, S, V, K, Tmax, T...
    %also, create independent functions for each method call...
    %stand alone functions in separate files

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
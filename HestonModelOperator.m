function [A,B] = HestonModelOperator(NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho)
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
    %the code above is the generator of A and B -- these two will be passed
    %as parameters to HestonMatvec2 that will be called LowRankMatVec

    A(:,:,1)=A1KX;
    A(:,:,2)=A2KX;
    A(:,:,3)=A3KX;
    A(:,:,4)=A4KX;
    A(:,:,5)=A5KX;

    B(:,:,1)=A1KY;
    B(:,:,2)=A2KY;
    B(:,:,3)=A3KY;
    B(:,:,4)=A4KY;
    B(:,:,5)=A5KY;
end

function [BX,BY] = HestonMatVecBoundaries(NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, K, Tmax, t, T)

    SMatrix = diag(S);
    S2Matrix = diag(S.^2);
    VMatrix = diag(V);

    b1x = zeros(NS,1);
    b1x(NS)=(1/ds);
    b1y = ones(NV,1);

    b2x = zeros(NS,1);
    b2x(NS)=(1/ds);
    b2y = ones(NV,1);
    
    discountedPayoff = max((S - K * exp(-r * (Tmax - T(t)))), 0);

    b3x = [discountedPayoff', S'];%In this case rank 2 because of the shape of the condition
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

    BX = [b1x,b2x,b3x,b4x,b5x];
    BY = [b1y,b2y,b3y,b4y,b5y];

end
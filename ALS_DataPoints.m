function [ClosedPrice, price, error, time] = ALS_DataPoints(nS, nV, nT, S0, V0, cnIterations)

    % Illustration of pricing using uniform and non-uniform grids

    
    % restoredefaultpath
    % rehash toolboxcache
    
    % Strike price, risk free rate, dividend yield, and maturity
    K = 100;
    r = 0.02;
    q = 0.05;
    Mat = 0.15;
    
    % Heston parameters
    kappa =  1.5;
    theta =  0.04;
    sigma =  0.3;
    rho   = -0.9;
    v0    =  0.05;
    lambda = 0;
    params = [kappa theta sigma v0 rho lambda];
    
    % Minimum and maximum values for the Stock Price, Volatility, and Maturity
    Smin = 0;  Smax = 2*K;
    Vmin = 0;  Vmax = 0.5;
    Tmin = 0;  Tmax = Mat;
    
    % Number of grid points for the stock, volatility, and maturity
    % nS = 79;        % Stock price
    % nV = 39;        % Volatility
    
    % nS = 29;        % Stock price -- 99 --29
    % nV = 19;        % Volatility  -- 99 --19
    % nT = 100;      % Maturity    -- 1000 
    
    % The maturity time increment and grid
    dt = (Tmax-Tmin)/nT;
    T = [0:nT].*dt;
    
    %% Pricing Using a Uniform Grid
    % Increment for Stock Price and volatility
    ds = (Smax-Smin)/nS;
    dv = (Vmax-Vmin)/nV;
    
    % The stock price and volatility grids
    S = [0:nS].*ds;
    V = [0:nV].*dv;
    
    Sm = [1:(nS-1)].*ds;
    Vm = [1:(nV-1)].*dv;
    
    
    epsilon = 0.00001;
    
    
    iterations = 10;
    restart = 3;
    
    % 
    % tic
    % UvHEClassicCNGMRS = HestonExplicitClassicCNRC01(params,K,r,q,Sm,Vm,T,1,iterations,restart);
    % timeElapsed = toc;
    
    
    
    % GMRES_Result = HestonExplicitClassicCN_GMRS(params,K,r,q,S,V,T);
    
    %P = HestonExplicitClassicCN_GMRS(params,K,r,q,S,V,T);
    % tic
    % UvHEClassicCNXYdev = HestonExplicitClassicCNXYRC04(params,K,r,q,Sm,Vm,T,2, iterations, restart);
    % timeElapsedXY = toc;
    
    tic
    UvHEClassicCNXYLeanDev = HestonExplicitClassicCNALSDev10(params,K,r,q,Sm,Vm,T,2,iterations, restart, cnIterations);
    %UvHEClassicCNXYLeanDev = HestonExplicitClassicCNXYRC04(params,K,r,q,Sm,Vm,T,2, iterations, restart);
    timeElapsedLean = toc;
    
    %Suliko = HestonExplicitClassicCNXYCOMP(params,K,r,q,S,V,T);
    
    % HestonExplicitClassicCN_GMRS(params,K,r,q,S,V,T)
    % test = HestonExplicitClassicCN_GMRS(params,K,r,q,Sm,Vm,T);
    
    % Obtain the price by 2-D interpolation
    % S0 = 101.52;
    % V0 = 0.05412;
    %UniformPrice = interp2(V,S,U,V0,S0);
    
    % UvHEClassicCNXYdevPrice = interp2(Vm,Sm,UvHEClassicCNXYdev,V0,S0);
    UvHEClassicCNXYLeanDevPrice = interp2(Vm,Sm,UvHEClassicCNXYLeanDev,V0,S0);
    
    %% Pricing Using a Non-Uniform Grid
    % The stock price grid
    c = K/5;
    dz = 1/nS*(asinh((Smax-K)/c) - asinh(-K/c));
    for i=1:nS+1;
	    z(i) = asinh(-K/c) + (i-1)*dz;
	    S(i) = K + c*sinh(z(i));
    end
    
    % The volatility grid
    d = Vmax/500;
    dn = asinh(Vmax/d)/nV;
    for j=1:nV+1
	    n(j) = (j-1)*dn;
	    V(j) = d*sinh(n(j));
    end
    
    % Solve the PDE
    %U = HestonExplicitPDENonUniformGrid(params,K,r,q,S,V,T);
    
    % Obtain the price by 2-D interpolation
    % NonUniformPrice = interp2(V,S,U,V0,S0);
    
    
    %% Closed form Price and errors
    trap = 1;
    PutCall = 'C';
    [x w] = GenerateGaussLaguerre(32);
    ClosedPrice = HestonPriceGaussLaguerre(PutCall,S0,K,Mat,r,q,kappa,theta,sigma,lambda,V0,rho,trap,x,w);
    
    % UvHEClassicCNXYdevError = UvHEClassicCNXYdevPrice - ClosedPrice;
    UvHEClassicCNXYLeanDevError = UvHEClassicCNXYLeanDevPrice - ClosedPrice;
    
    
    % %% Output the results
    % %clc;
    % fprintf('Stock price grid size  %5.0f\n', nS+1)
    % fprintf('Volatility grid size   %5.0f\n', nV+1)
    % fprintf('Number of time steps   %5.0f\n', nT)
    % fprintf('-----------------------------------------------------------------------------------\n')
    % fprintf('Method                                          Price    DollarError Execution time\n')
    % fprintf('-----------------------------------------------------------------------------------\n')
    % fprintf('Closed Form                                %10.4f              \n', ClosedPrice)
    % fprintf('Heston Classic CN GMRS Low Rank Super Dima %10.4f        %5.2f       %10.4f\n', UvHEClassicCNXYdevPrice,UvHEClassicCNXYdevError, timeElapsedXY)
    % %fprintf('Heston Classic CN GMRS Low Rank Lean       %10.4f        %5.2f       %10.4f\n\n', UvHEClassicCNXYLeanDevPrice,UvHEClassicCNXYLeanDevError, timeElapsedLean)
    % fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %5.2f       %10.4f\n\n', UvHEClassicCNXYLeanDevPrice,UvHEClassicCNXYLeanDevError, timeElapsedLean)
    % %fprintf('Dmitry is the BEST                                                         \n')
    % fprintf('-----------------------------------------------------------------------------------\n')
    
    price = UvHEClassicCNXYLeanDevPrice;
    error = UvHEClassicCNXYLeanDevError;
    time = timeElapsedLean;

end

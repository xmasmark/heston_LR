% Illustration of pricing using uniform and non-uniform grids
clc; clear;


% HestonExplicitClassicCN_GMRS(params,K,r,q,S,V,T)
% test = HestonExplicitClassicCN_GMRS(params,K,r,q,Sm,Vm,T);

% Obtain the price by 2-D interpolation
S0 = 101.52;
V0 = 0.05412;
%UniformPrice = interp2(V,S,U,V0,S0);

%[price, error, time] = ALS_DataPoints(nS, nV, nT, S0, V0, cnIterations);
[price, error, time] = ALS_DataPoints(29, 19, 100, S0, V0, 100);

%% Output the results
%clc;
% fprintf('Stock price grid size  %5.0f\n', nS+1)
% fprintf('Volatility grid size   %5.0f\n', nV+1)
% fprintf('Number of time steps   %5.0f\n', nT)
fprintf('-----------------------------------------------------------------------------------\n')
fprintf('Method                                          Price    DollarError Execution time\n')
fprintf('-----------------------------------------------------------------------------------\n')
% fprintf('Closed Form                                %10.4f              \n', ClosedPrice)
% fprintf('Uniform Grid Book %10.4f    %5.2f\n', UniformPrice,UError)
% fprintf('Non-Uniform Grid  %10.4f    %5.2f\n', NonUniformPrice,NError)
% fprintf('Marco Heston Algo %10.4f    %5.2f\n', UniformPriceLC,MError)
% fprintf('Marco Heston LR   %10.4f    %5.2f\n', UniformPriceLowRank,LRError)
% fprintf('Marco Heston UVT                 %10.4f    %5.2f\n', UniformPriceTilda,UVTildaError)
% fprintf('Marco Heston UVT Plus            %10.4f    %5.2f\n', UniformPriceTildaPlus,UVTildaPlusError)
% fprintf('Marco Heston UVT GMRES           %10.4f    %5.2f\n', UniformPriceTildaPlusGMRES,UVTildaPlusErrorGMRES)
% fprintf('Marco Heston UVT3 %10.4f    %5.2f\n', UniformPriceTilda3,UVTildaError3)
% fprintf('Marco Heston Euler               %10.4f    %5.2f\n', UvHestonEulerPrice,UvHestonEulerError)
%fprintf('Marco Heston CN        %10.4f    %5.2f\n', UvHestonCNPrice,UvHestonCNError)
% fprintf('Heston Classic Euler             %10.4f    %5.2f\n', UniformPriceHEClassicEuler,UvHEClassicEulerError)
% fprintf('Heston Classic Euler Operator              %10.4f          %5.2f\n', UniformPriceHEClassicEulerIntegrated,UvHEClassicEulerIntegratedError)
% fprintf('Heston Classic CN                %10.4f    %5.2f\n', UniformPriceHEClassicCN,UvHEClassicCNError)
% fprintf('Heston Classic CN GMRS vectorised          %10.4f        %5.2f       %10.4f\n', UniformPriceHEClassicCNGMRS,UvHEClassicCNGMRSError, timeElapsed)
%fprintf('Heston Classic CN Backslash vectorised     %10.4f        %5.2f       %10.4f\n', UniformPriceHEClassicCNGMRS,UvHEClassicCNGMRSError, timeElapsed)
% fprintf('Heston Classic CN GMRS Low Rank Super Dima %10.4f        %5.2f       %10.4f\n', UvHEClassicCNXYdevPrice,UvHEClassicCNXYdevError, timeElapsedXY)
%fprintf('Heston Classic CN GMRS Low Rank Lean       %10.4f        %5.2f       %10.4f\n\n', UvHEClassicCNXYLeanDevPrice,UvHEClassicCNXYLeanDevError, timeElapsedLean)
fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %5.2f       %10.4f\n\n', price,error, time)
%fprintf('Dmitry is the BEST                                                         \n')
fprintf('-----------------------------------------------------------------------------------\n')

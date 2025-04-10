% Illustration of pricing using uniform and non-uniform grids
clc; clear;


% HestonExplicitClassicCN_GMRS(params,K,r,q,S,V,T)
% test = HestonExplicitClassicCN_GMRS(params,K,r,q,Sm,Vm,T);

% Obtain the price by 2-D interpolation
S0 = 101.52;
V0 = 0.05412;
%UniformPrice = interp2(V,S,U,V0,S0);

%[price, error, time] = ALS_DataPoints(nS, nV, nT, S0, V0, cnIterations);
nt = 1000;
% % % % % % % [ClosedPrice10, price10, error10, time10] = ALS_DataPoints(29, 19, nt, S0, V0, 10);
% % % % % % % [ClosedPrice20, price20, error20, time20] = ALS_DataPoints(29, 19, nt, S0, V0, 20);
% % % % % % % [ClosedPrice30, price30, error30, time30] = ALS_DataPoints(29, 19, nt, S0, V0, 30);
% % % % % % % [ClosedPrice40, price40, error40, time40] = ALS_DataPoints(29, 19, nt, S0, V0, 40);
% % % % % % % [ClosedPrice50, price50, error50, time50] = ALS_DataPoints(29, 19, nt, S0, V0, 50);
% % % % % % % [ClosedPrice100, price100, error100, time100] = ALS_DataPoints(29, 19, nt, S0, V0, 100);
% % % % % % % % [ClosedPrice200, price200, error200, time200] = ALS_DataPoints(29, 19, 100, S0, V0, 200);
% % % % % % % % [ClosedPrice300, price300, error300, time300] = ALS_DataPoints(29, 19, 100, S0, V0, 300);
% % % % % % % 
% % % % % % % %% Output the results
% % % % % % % fprintf('-------------------------------------------------------------------------------------------------------------------------------------\n')
% % % % % % % fprintf('Method                                          Price      Closed Price DollarError  Execution time       NS       NV    CNIterations\n')
% % % % % % % fprintf('-------------------------------------------------------------------------------------------------------------------------------------\n')
% % % % % % % fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price10, ClosedPrice10,error10, time10, 30, 20, 10)
% % % % % % % fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price20, ClosedPrice20,error20, time20, 30, 20, 20)
% % % % % % % fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price30, ClosedPrice20,error30, time30, 30, 20, 30)
% % % % % % % fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price40, ClosedPrice20,error40, time40, 30, 20, 40)
% % % % % % % fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price50, ClosedPrice20,error50, time50, 30, 20, 50)
% % % % % % % fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price100, ClosedPrice100,error100, time100, 30, 20, 100)
% % % % % % % % fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price200, ClosedPrice200,error200, time200, 30, 20, 200)
% % % % % % % % fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price300, ClosedPrice300,error300, time300, 30, 20, 300)
% % % % % % % %fprintf('Dmitry is the BEST                                                         \n')
% % % % % % % fprintf('-------------------------------------------------------------------------------------------------------------------------------------\n')
% % % % % % % 


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

             % HestonVanillaClosedForm(S,  K,   T, r, kappa, theta, sigma, rho, v0, type)
closePriceRL = HestonVanillaClosedForm(S0, K, Mat, r, kappa, theta, sigma, rho, v0, 'call');

%[price, error, time] = ALS_DataPoints(nS, nV, nT, S0, V0, cnIterations);
nt = 200;
cnIterations = 50;
ns=99;
nv=99;
[ClosedPrice10, price10, error10, time10] = ALS_DataPoints(ns, nv, nt, S0, V0, cnIterations);
% [ClosedPrice20, price20, error20, time20] = ALS_DataPoints(39, 29, nt, S0, V0, cnIterations);
% [ClosedPrice30, price30, error30, time30] = ALS_DataPoints(49, 39, nt, S0, V0, cnIterations);
% [ClosedPrice40, price40, error40, time40] = ALS_DataPoints(59, 49, nt, S0, V0, cnIterations);
% [ClosedPrice50, price50, error50, time50] = ALS_DataPoints(69, 59, nt, S0, V0, cnIterations);
% [ClosedPrice100, price100, error100, time100] = ALS_DataPoints(79, 69, nt, S0, V0, cnIterations);
% [ClosedPrice200, price200, error200, time200] = ALS_DataPoints(29, 19, 100, S0, V0, 200);
% [ClosedPrice300, price300, error300, time300] = ALS_DataPoints(29, 19, 100, S0, V0, 300);

%% Output the results
fprintf('-------------------------------------------------------------------------------------------------------------------------------------\n')
fprintf('Method                                          Price      Closed Price DollarError  Execution time       NS       NV    CNIterations\n')
fprintf('-------------------------------------------------------------------------------------------------------------------------------------\n')
fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price10, ClosedPrice10,error10, time10, ns, nv, cnIterations)
% fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price20, ClosedPrice20,error20, time20, 40, 30, 100)
% fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price30, ClosedPrice20,error30, time30, 50, 40, 100)
% fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price40, ClosedPrice20,error40, time40, 60, 50, 100)
% fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price50, ClosedPrice20,error50, time50, 70, 60, 100)
% fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price100, ClosedPrice100,error100, time100, 80, 70, 100)
% fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price200, ClosedPrice200,error200, time200, 30, 20, 200)
% fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price300, ClosedPrice300,error300, time300, 30, 20, 300)
%fprintf('Dmitry is the BEST                                                         \n')
fprintf('-------------------------------------------------------------------------------------------------------------------------------------\n')

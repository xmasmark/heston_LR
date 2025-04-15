% Illustration of pricing using uniform and non-uniform grids
clc; clear;


% HestonExplicitClassicCN_GMRS(params,K,r,q,S,V,T)
% test = HestonExplicitClassicCN_GMRS(params,K,r,q,Sm,Vm,T);

% Obtain the price by 2-D interpolation
S0 = 101.52;
V0 = 0.05412;
%UniformPrice = interp2(V,S,U,V0,S0);

%[price, error, time] = ALS_DataPoints(nS, nV, nT, S0, V0, cnIterations);
% nt = 1000;
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

%[price, error, time] = ALS_DataPoints(nS, nV, nT, S0, V0, cnIterations);
nt = 200;
cnIterations = 50;


%[ClosedPrice, price, error, time] = GMRES_Low_Rank_DataPoints(nS, nV, nT, S0, V0, cnIterations)

% [ClosedPrice10, price10, error10, time10] = GMRES_Low_Rank_DataPoints(69, 59, nt, S0, V0, cnIterations);
% [ClosedPrice20, price20, error20, time20] = GMRES_Low_Rank_DataPoints(99, 89, nt, S0, V0, cnIterations);
% [ClosedPrice30, price30, error30, time30] = GMRES_Low_Rank_DataPoints(129, 119, nt, S0, V0, cnIterations);
% [ClosedPrice40, price40, error40, time40] = GMRES_Low_Rank_DataPoints(159, 149, nt, S0, V0, cnIterations);
%[ClosedPrice50, price50, error50, time50] = GMRES_Low_Rank_DataPoints(49, 49, nt, S0, V0, cnIterations);
%[ClosedPrice100, price100, error100, time100] = GMRES_Low_Rank_DataPoints(299, 299, nt, S0, V0, cnIterations);

standardNS = 799;
standardNV = 49;
% nt200 = 200;
% [ClosedPrice200, price200, error200, time200] = GMRES_Low_Rank_DataPoints(standardNS, standardNV, nt200, S0, V0, cnIterations);
% nt300 = 300;
% [ClosedPrice300, price300, error300, time300] = GMRES_Low_Rank_DataPoints(standardNS, standardNV, nt300, S0, V0, cnIterations);
% nt400 = 400;
% [ClosedPrice400, price400, error400, time400] = GMRES_Low_Rank_DataPoints(standardNS, standardNV, nt400, S0, V0, cnIterations);
nt500 = 1000;
[ClosedPrice500, price500, error500, time500] = GMRES_Low_Rank_DataPoints(standardNS, standardNV, nt500, S0, V0, cnIterations);
% [ClosedPrice300, price300, error300, time300] = ALS_DataPoints(29, 19, 100, S0, V0, 300);

%% Output the results
fprintf('-------------------------------------------------------------------------------------------------------------------------------------\n')
fprintf('Method                                          Price      Closed Price DollarError  Execution time       NS       NV    nt\n')
fprintf('-------------------------------------------------------------------------------------------------------------------------------------\n')
% fprintf('Heston Classic GMRES Low Rank Lean         %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price10, ClosedPrice10,error10, time10, 70, 60, 100)
% fprintf('Heston Classic GMRES Low Rank Lean         %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price20, ClosedPrice20,error20, time20, 100, 90, 100)
% fprintf('Heston Classic GMRES Low Rank Lean         %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price30, ClosedPrice20,error30, time30, 130, 120, 100)
% fprintf('Heston Classic GMRES Low Rank Lean         %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price40, ClosedPrice20,error40, time40, 160, 150, 100)
%fprintf('Heston Classic GMRES Low Rank Lean         %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price50, ClosedPrice50,error50, time50, 500, 500, 100)
%fprintf('Heston Classic GMRES Low Rank Lean         %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price100, ClosedPrice100,error100, time100, 300, 300, 100)

% fprintf('Heston Classic GMRES Low Rank Lean         %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price200, ClosedPrice200, error200, time200, 300, 300, nt200)
% 
% fprintf('Heston Classic GMRES Low Rank Lean         %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price300, ClosedPrice300, error300, time300, 300, 300, nt300)
% 
% fprintf('Heston Classic GMRES Low Rank Lean         %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price400, ClosedPrice400, error400, time400, 300, 300, nt400)

fprintf('Heston Classic GMRES Low Rank Lean         %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price500, ClosedPrice500, error500, time500, standardNS+1, standardNV+1, nt500)
% fprintf('Heston Classic CN ALS  Low Rank Lean       %10.4f        %10.4f      %5.2f       %10.4f       %d       %d             %d\n',price300, ClosedPrice300,error300, time300, 30, 20, 300)
%fprintf('Dmitry is the BEST                                                         \n')
fprintf('-------------------------------------------------------------------------------------------------------------------------------------\n')

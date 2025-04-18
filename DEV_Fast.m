% Illustration of pricing using uniform and non-uniform grids
clc; clear;

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

nS = 29;        % Stock price
nV = 19;        % Volatility
nT = 200;      % Maturity

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

% Solve the PDE
%U = HestonExplicitPDE(params,K,r,q,S,V,T);

% Un = HestonExplicitPDELinearComplexity(params,K,r,q,S,V,T);
% 
% Uvector = HestonExplicitWithVectors(params,K,r,q,S,V,T);
% 
% UvLR = HestonExplicitLowRank(params,K,r,q,S,V,T);

epsilon = 0.00001;

% numIterations = 1;
% 
% % Preallocate array to store execution times
% execution_times = zeros(1, numIterations);
% 
% % Loop to call the function and measure execution time
% for i = 1:numIterations
%     tic; % Start timing
%     UvLR = HestonExplicitLowRank(params,K,r,q,S,V,T);
%     execution_times(i) = toc; % Stop timing and store the elapsed time
% end
% 
% % Display the execution times (optional)
% %disp(execution_times);
% 
% % Optionally, you can calculate and display some statistics
% avg_time = min(execution_times);
% disp(['Minimum execution time of HestonExplicit PDE - Marco Implementation: ', num2str(avg_time), ' seconds']);
% 
% execution_times = zeros(1, numIterations);
% 
% for i = 1:numIterations
%     tic; % Start timing
%     %profile on
%     UvTilda3 = HestonExplicitLowRankRTildaRC3(params,K,r,q,S,V,T,epsilon);
%     %profile viewer
%     execution_times(i) = toc; % Stop timing and store the elapsed time
% end

%UCN = HestonExplicitCrankNicholson(params,K,r,q,S,V,T,epsilon);

%UvHestonEuler = HestonExplicitEuler(params,K,r,q,S,V,T,epsilon);

% % Display the execution times (optional)
% %disp(execution_times);
% 
% % Optionally, you can calculate and display some statistics
% avg_time = min(execution_times);
% disp(['Minimum execution time of HestonExplicit LowRank R TildaRC3: ', num2str(avg_time), ' seconds']);
% 
%UvTilda = HestonExplicitLowRankRTildaRC3(params,K,r,q,S,V,T,epsilon);

%UvTildaPlus = HestonExplicitLowRankRTildaPlus(params,K,r,q,S,V,T,epsilon);
%UvTildaPlusGMRES = HestonExplicitLowRankRTildaPlusGMRES(params,K,r,q,S,V,T,epsilon);

%HestonExplicitClassicXYEuler(params,K,r,q,S,V,T,mode)

%HestonExplicitClassic(params,K,r,q,S,V,T)
%UvHEClassicEuler = HestonExplicitClassic(params,K,r,q,Sm,Vm,T);
UvHEClassicEulerIntegrated = HestonExplicitClassicXYEuler(params,K,r,q,Sm,Vm,T,0);

%HestonExplicitClassicCNRC01
%UvHEClassicCN = HestonExplicitClassicCN(params,K,r,q,Sm,Vm,T);
%UvHEClassicCN = HestonExplicitClassicCNRC01(params,K,r,q,Sm,Vm,T,1);

%UvHEClassicCNGMRS = HestonExplicitClassicCN_GMRS(params,K,r,q,Sm,Vm,T);

iterations = 10;
restart = 3;

UvHEClassicCNGMRS = HestonExplicitClassicCNRC01(params,K,r,q,Sm,Vm,T,2, iterations, restart);

% GMRES_Result = HestonExplicitClassicCN_GMRS(params,K,r,q,S,V,T);

%P = HestonExplicitClassicCN_GMRS(params,K,r,q,S,V,T);

UvHEClassicCNXYdev = HestonExplicitClassicCNXYRC02(params,K,r,q,Sm,Vm,T,2, iterations, restart);
                     %HestonExplicitClassicCNALSDev01(params,K,r,q,S,V,T, mode, iterations, restart)
UvHEClassicCNALSdev = HestonExplicitClassicCNALSDev01(params,K,r,q,Sm,Vm,T,2,iterations, restart);

%Suliko = HestonExplicitClassicCNXYCOMP(params,K,r,q,S,V,T);

% HestonExplicitClassicCN_GMRS(params,K,r,q,S,V,T)
% test = HestonExplicitClassicCN_GMRS(params,K,r,q,Sm,Vm,T);

% Obtain the price by 2-D interpolation
S0 = 101.52;
V0 = 0.05412;
%UniformPrice = interp2(V,S,U,V0,S0);

% % % S0 = 101.52;
% % % V0 = 0.05412;
% % % %UniformPriceLC = interp2(V,S,Un,V0,S0);
% % % 
% % % S0 = 101.52;
% % % V0 = 0.05412;
% % % %UniformPriceLowRank = interp2(V,S,UvLR,V0,S0);

%UniformPriceTilda = interp2(V,S,UvTilda,V0,S0);
%UniformPriceTilda3 = interp2(V,S,UvTilda3,V0,S0);
%UvHestonEulerPrice = interp2(V,S,UvHestonEuler,V0,S0);
%UvHestonCNPrice = interp2(V,S,UCN,V0,S0);

%UniformPriceTildaPlus = interp2(V,S,UvTildaPlus,V0,S0);
%UniformPriceTildaPlusGMRES = interp2(V,S,UvTildaPlusGMRES,V0,S0);

%UniformPriceHEClassicEuler = interp2(Vm,Sm,UvHEClassicEuler,V0,S0);

UniformPriceHEClassicEulerIntegrated = interp2(Vm,Sm,UvHEClassicEulerIntegrated,V0,S0);

%UniformPriceHEClassicCN = interp2(Vm,Sm,UvHEClassicCN,V0,S0);

UniformPriceHEClassicCNGMRS = interp2(Vm,Sm,UvHEClassicCNGMRS,V0,S0);

UvHEClassicCNXYdevPrice = interp2(Vm,Sm,UvHEClassicCNXYdev,V0,S0);

UvHEClassicCNALSdevPrice = interp2(Vm,Sm,UvHEClassicCNALSdev,V0,S0);

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
U = HestonExplicitPDENonUniformGrid(params,K,r,q,S,V,T);

% Obtain the price by 2-D interpolation
NonUniformPrice = interp2(V,S,U,V0,S0);


%% Closed form Price and errors
trap = 1;
PutCall = 'C';
[x w] = GenerateGaussLaguerre(32);
ClosedPrice = HestonPriceGaussLaguerre(PutCall,S0,K,Mat,r,q,kappa,theta,sigma,lambda,V0,rho,trap,x,w);
%UError = UniformPrice - ClosedPrice;
%NError = NonUniformPrice  - ClosedPrice;
%MError = UniformPriceLC - ClosedPrice;
%LRError = UniformPriceLowRank - ClosedPrice;
% UVTildaError = UniformPriceTilda - ClosedPrice;
% UVTildaPlusError = UniformPriceTildaPlus - ClosedPrice;
% UVTildaPlusErrorGMRES = UniformPriceTildaPlusGMRES - ClosedPrice;
%UVTildaError3 = UniformPriceTilda3 - ClosedPrice;
% UvHestonEulerError = UvHestonEulerPrice - ClosedPrice;
% UvHestonCNError = UvHestonCNPrice - ClosedPrice;
% UvHEClassicEulerError = UniformPriceHEClassicEuler - ClosedPrice;
UvHEClassicEulerIntegratedError = UniformPriceHEClassicEulerIntegrated - ClosedPrice;
% UvHEClassicCNError = UniformPriceHEClassicCN - ClosedPrice;
UvHEClassicCNGMRSError = UniformPriceHEClassicCNGMRS - ClosedPrice;
UvHEClassicCNXYdevError = UvHEClassicCNXYdevPrice - ClosedPrice;
UvHEClassicCNALSdevError = UvHEClassicCNALSdevPrice - ClosedPrice;


%% Output the results
%clc;
fprintf('Stock price grid size  %5.0f\n', nS+1)
fprintf('Volatility grid size   %5.0f\n', nV+1)
fprintf('Number of time steps   %5.0f\n', nT)
fprintf('---------------------------------------------------------------------------\n')
fprintf('Method                                          Price    DollarError       \n')
fprintf('---------------------------------------------------------------------------\n')
fprintf('Closed Form                                %10.4f              \n', ClosedPrice)
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
fprintf('Heston Classic Euler Operator              %10.4f          %5.2f\n', UniformPriceHEClassicEulerIntegrated,UvHEClassicEulerIntegratedError)
% fprintf('Heston Classic CN                %10.4f    %5.2f\n', UniformPriceHEClassicCN,UvHEClassicCNError)
fprintf('Heston Classic CN GMRS vectorised          %10.4f          %5.2f\n', UniformPriceHEClassicCNGMRS,UvHEClassicCNGMRSError)
fprintf('Heston Classic CN GMRS Low Rank            %10.4f          %5.2f\n', UvHEClassicCNXYdevPrice,UvHEClassicCNXYdevError)
fprintf('Heston Classic CN ALS  Low Rank            %10.4f          %5.2f\n', UvHEClassicCNALSdevPrice,UvHEClassicCNALSdevError)
fprintf('---------------------------------------------------------------------------\n')

function U = HestonExplicitPDELinearComplexity(params,K,r,q,S,V,T)

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
epsilon = 0.01;

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
for s=1:NS
	for v=1:NV
		U(s,v) = max(S(s) - K, 0);
	end
end

matrix_rank=[];
matrix_rank_home_made=[];

for t=1:NT-1
	% Boundary condition for Smin and Smax
	for v=1:NV-1
		U(1,v) = 0;
		U(NS,v) = max(0, Smax - K); % Galiotos uses U(NS-1,v) + ds;
	end
	% Boundary condition for Vmax
	for s=1:NS
		U(s,NV) = max(0, S(s) - K); % Galiotos uses U(s,NV-1);
	end
	% Update the temporary grid u(s,t) with the boundary conditions
	u = U;
	% Boundary condition for Vmin.
	% Previous time step values are in the temporary grid u(s,t)
	for s=2:NS-1
		DerV = (u(s,2) - u(s,1)) / dv;                	% PDE Points on the middle of the grid (non boundary)
		DerS = (u(s+1,1) - u(s-1,1))/2/ds;              % Central difference for dU/dS
		U(s,1) = u(s,1)*(1 - r*dt - kappa*theta*dt/dv)...
			   + dt*(0.5*(r-q)*s*u(s+1,1) - u(s-1,1)) ...
			   + kappa*theta*dt/dv*u(s,2);
	end
	% Update the temporary grid u(s,t) with the boundary conditions
	u = U;
	% Interior points of the grid (non boundary).
	% Previous time step values are in the temporary grid u(s,t)

	% for s=2:NS-1
	% 	for v=2:NV-1
		% 	A = (1 - dt*(s-1)^2*(v-1)*dv - sigma^2*(v-1)*dt/dv - r*dt);
		% 	B = (1/2*dt*(s-1)^2*(v-1)*dv - 1/2*dt*(r-q)*(s-1));
		% 	C = (1/2*dt*(s-1)^2*(v-1)*dv + 1/2*dt*(r-q)*(s-1));
		% 	D = (1/2*dt*sigma^2*(v-1)/dv - 1/2*dt*kappa*(theta-(v-1)*dv)/dv);
		% 	E = (1/2*dt*sigma^2*(v-1)/dv + 1/2*dt*kappa*(theta-(v-1)*dv)/dv);
		% 	F = 1/4*dt*sigma*(s-1)*(v-1);
		% 	U(s,v) = A*u(s,v) + B*u(s-1,v) + C*u(s+1,v)...   % The PDE
			% 	   + D*u(s,v-1) + E*u(s,v+1)...
			% 	   + F*(u(s+1,v+1)+u(s-1,v-1)-u(s-1,v+1)-u(s+1,v-1));
	% 	end
    % end
    % remember U indexation: U=U(s,v)
    % check page 317 in the book - this is explained
	for s=2:NS-1
		for v=2:NV-1
            %Here I will write the various components of the PDE

            % derS=(u(s+1,v)-u(s,v))/ds;
            % derS2=(u(s+1,v)-2*u(s,v)+u(s-1,v))/(ds^2);
            % derV=(u(s,v+1)-u(s,v))/dv;
            % derV2=(u(s,v+1)-2*u(s,v)+u(s,v-1))/(dv^2);
            % derVS=(u(s+1,v+1)+u(s-1,v-1)-u(s-1,v+1)-u(s+1,v-1))/(4*dv*ds);

            %use central approximation rather than forward

            derS=(u(s+1,v)-u(s-1,v))/(2*ds);
            derS2=(u(s+1,v)-2*u(s,v)+u(s-1,v))/(ds^2);
            derV=(u(s,v+1)-u(s,v-1))/(2*dv);
            derV2=(u(s,v+1)-2*u(s,v)+u(s,v-1))/(dv^2);
            derVS=(u(s+1,v+1)+u(s-1,v-1)-u(s-1,v+1)-u(s+1,v-1))/(4*dv*ds);

            pa = u(s,v);
            pb = u(s,v+1);
            pc = u(s,v-1);

            pd = u(s+1,v);
            pe = u(s+1,v+1);
            pf = u(s+1,v-1);

            pg = u(s-1,v);
            ph = u(s-1,v+1);
            pi = u(s-1,v-1);

            %U(s,v)=0.5*v*s^2*derS2+rho*sigma*v*s*derVS+0.5*sigma^2*v*derV2-r*u(s,v)+(r-q)*s*derS+kappa*(theta-v)*derV;
            %The following is formula 10.19 at page 316 of the book
            %U(s,v)=u(s,v)+dt*u(s,v)*((0.5*v*s^2*derS2)+(sigma*v*s*derVS)+(0.5*sigma^2*v*derV2)+((r-q)*s*derS)+kappa*(theta-v)*derV-r);
            %U(s,v)=u(s,v)+dt*u(s,v)*((0.5*v*s^2*derS2)+(0.5*sigma^2*v*derV2)+(sigma*v*s*derVS)+((r-q)*s*derS)+(kappa*(theta-v)*derV)-r);

            sValue =S(s);
            vValue=V(v);

            old_value = u(s,v);

            first = (0.5*vValue*sValue^2*derS2);
            second = (0.5*sigma^2*vValue*derV2);
            third = (sigma*vValue*sValue*derVS);
            fourth = ((r-q)*sValue*derS);
            fifth = (kappa*(theta-vValue)*derV)-r*u(s,v);

            U(s,v)=u(s,v)+dt*(first+second+third+fourth+fifth);

            current = U(s,v);
            if(isnan(current))
                aa= current;
            end  
		end
    end


    % builtin_rank = rank(U,epsilon);
    % 
    % [uc,sc,vs]=svd(U);
    % 
    % sigSum = 0;
    % 
    % %the look should go to the minimum of the two dimensions
    % %so a step is missing and it should be added when moving
    % %this to a separate function
    % for i = 1: size(sc,2)
    %     sigSum = sigSum + power(sc(i,i),2);
    % end
    % 
    % threshold = power(epsilon,2)*sigSum;
    % 
    % rank_home_made = 0;
    % 
    % errorSquare = 0;
    % 
    % for i=size(sc,2):-1:1
    %     errorSquare = errorSquare + power(sc(i,i),2);
    %     if errorSquare > threshold
    %         rank_home_made = i + 1;
    %         break;
    %     end
    % end
    % 
    % % if t==2900
    % %     [X,Y]=meshgrid(1:NV,1:NS);
    % %     surf(X,Y,U)
    % % end
    % 
    % 
    % matrix_rank_home_made = [t,s,v,rank(U,0.01),rank_home_made;matrix_rank_home_made];
    % 
    % %creation of approximated matrix
    % %[uc,sc,vs]=svd(U);
    % 
    % if rank_home_made > 0 && t == 1500
    % 
    %     ucshort = uc;
    %     ucshort(:,rank_home_made+1:size(ucshort,2))=0;
    %     vsshort = vs;
    %     vsshort(rank_home_made+1:size(vsshort,1),:)=0;
    %     vsshort(rank_home_made+1:size(vsshort,1),:)=0;
    % 
    %     scshort = sc;
    % 
    %     for i = rank_home_made+1 : size(sc,2)
    %         scshort(i,i)=0;
    %     end
    % 
    %     Uapproximated=ucshort*scshort*vsshort;
    % 
    %     DeltaMatrix = U-Uapproximated;
    %     n = norm(DeltaMatrix,"fro");
    % 
    % end
    % 

    %plot at time zero, time 1500, time 3000

end

%writematrix(matrix_rank, 'G:\My Drive\BO\PhD\Heston\RankAnalysis\matrix_rank.csv')
%writematrix(matrix_rank_home_made, 'G:\My Drive\BO\PhD\Heston\RankAnalysis\matrix_rank_home_made.csv')

%Task: rewrite the code for 1 step only (t=0->t=1)
%In the new code I don't use 2 dimensional U and u: it has to work with
%2 vectors that are 1 dimensional (1 row and 1 column).
%the start is U that is 2 vectors, 1 depending on S and one depending on V.
%Like it was done for the boundaries. The stuff that is currently written
%in the two dimensional case, needs to work on vectors rather than 2
%dimensional matrix. One step implementation will show the linear
%complexity leading to the same result of the quadratic complexity (current
%implementation).


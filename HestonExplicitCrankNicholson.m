function U = HestonExplicitCrankNicholson(params,K,r,q,S,V,T,epsilonError)
    % Finite differences for the Heston PDE for a European Call

    % Heston parameters
    kappa = params(1);
    theta = params(2);
    sigma = params(3);
    v0    = params(4);
    rho   = params(5);
    lambda = params(6);
    
    % Length of stock price, volatility, and maturity
    NS = length(S);
    NV = length(V);
    NT = length(T);
    Smin = S(1);  Smax = S(NS);
    Vmin = V(1);  Vmax = V(NV);
    Tmin = T(1);  Tmax = T(NT);
    
    % Grid spacing
    ds = (Smax-Smin)/(NS-1);
    dt = (Tmax-Tmin)/(NT-1);
    dv = (Vmax-Vmin)/(NV-1);

    % Initialize the solution matrix for a given time snapshot
    U = zeros(NS, NV);
    
    % Initial conditions: Call option payoff
    for i = 1:NS
        for j = 1:NV
            U(i,j) = max(S(i) - K, 0);
        end
    end

    u = zeros(NS, NV);

    % Crank-Nicholson coefficients and matrix setup
    for t = 1:NT
        % Boundary condition for Smin and Smax
        for v = 1:NV
            U(1,v) = 0;
            U(NS,v) = max(0, Smax - K);
        end
        % Boundary condition for Vmax
        for s = 1:NS
            U(s,NV) = max(0, S(s) - K);
        end
        % Update the temporary grid u(s,t) with the boundary conditions
        u = U;
        % Boundary condition for Vmin
        for s = 2:NS-1
            U(s,1) = max(u(s,1)*(1 - r*dt - kappa*theta*dt/dv) + dt*(0.5*(r-q)*s*u(s+1,1) - u(s-1,1)) + kappa*theta*dt/dv*u(s,2), 0);
        end

        % Construct the tridiagonal matrix for each j
        for j = 2:NV-1
            A = zeros(NS-2, NS-2);
            d = zeros(NS-2, 1);

            for i = 2:NS-1
                alpha   = 0.25 * dt * (sigma^2 * V(j) * (S(i)^2) - kappa * (theta - V(j)) * S(i));
                beta    = 0.25 * dt * (sigma^2 * V(j) * (S(i)^2) + kappa * (theta - V(j)) * S(i));
                gamma   = 0.25 * dt * (r * S(i) - rho * sigma * V(j) * S(i) * j);
                delta   = 0.25 * dt * (r * S(i) + rho * sigma * V(j) * S(i) * j);
                epsilon =  0.5 * dt * (r * V(j) - rho * sigma * V(j) * j);

                % Coefficients for the tridiagonal matrix
                a = -alpha - delta - epsilon;
                b = 1 + dt * r + 2 * epsilon;
                c = -beta - gamma - epsilon;

                % Add a small value to diagonal to avoid singularity
                b = b + 1e-10;

                % Fill the tridiagonal matrix A
                if i > 2
                    A(i-1, i-2) = a;
                end
                A(i-1, i-1) = b;  % This should be a scalar
                if i < NS-1
                    A(i-1, i) = c;
                end

                % Construct the right-hand side vector
                d(i-1) = alpha * U(i-1,j-1) + delta * U(i-1,j+1) + (1 - dt * r - 2 * epsilon) * U(i,j) + gamma * U(i+1,j-1) + beta * U(i+1,j+1);
            end

            % Solve the linear system using a more robust method
            u(2:NS-1, j) = A \ d;
        end

        % Ensure no negative values
        u = max(u, 0);
        
        % Update the solution matrix
        U = u;
    end
end



% % % function U = HestonExplicitCrankNicholson(params,K,r,q,S,V,T,epsilonError)
% % % 
% % %     % Finite differences for the Heston PDE for a European Call
% % %     % Uses even grid sizes
% % % 
% % %     % Heston parameters
% % %     kappa = params(1);
% % %     theta = params(2);
% % %     sigma = params(3);
% % %     v0    = params(4);
% % %     rho   = params(5);
% % %     lambda = params(6);
% % % 
% % %     % Length of stock price, volatility, and maturity
% % %     NS = length(S);
% % %     NV = length(V);
% % %     NT = length(T);
% % %     Smin = S(1);  Smax = S(NS);
% % %     Vmin = V(1);  Vmax = V(NV);
% % %     Tmin = T(1);  Tmax = T(NT);
% % % 
% % %     % Grid spacing
% % %     %dt = T / NT;
% % %     %ds = Smax / NS;
% % %     %dV = Vmax / NV;
% % %     ds = (Smax-Smin)/(NS-1);
% % %     dt = (Tmax-Tmin)/(NT-1);
% % %     dv = (Vmax-Vmin)/(NV-1);
% % % 
% % % 
% % %     % Grid spacing
% % %     %dt = T / NT;
% % %     dt = (Tmax-Tmin)/(NT-1);
% % % 
% % %     dS = Smax / NS;
% % %     dV = Vmax / NV;
% % % 
% % %     % Initialize the solution matrix for a given time snapshot
% % %     U = zeros(NS, NV);
% % % 
% % %     % Stock price and volatility vectors
% % %     %S = linspace(0, Smax, NS);
% % %     %V = linspace(0, Vmax, NV);
% % % 
% % %     % Boundary and initial conditions
% % %     % Example: Call option payoff
% % %     for i = 1:NS
% % %         for j = 1:NV
% % %             U(i,j) = max(S(i) - K, 0);
% % %         end
% % %     end
% % % 
% % %     u = zeros(NS, NV);
% % % 
% % %     % Crank-Nicolson coefficients and matrix setup
% % %     for t = 1:NT
% % % 
% % %         % % Enforce boundary conditions
% % %         % U_new(1, :) = 0; % S = 0
% % %         % U_new(NS, :) = Smax - K * exp(-r * (t * dt)); % S = Smax
% % %         % U_new(:, 1) = 0; % V = 0
% % %         % U_new(:, NV) = U(:, NV-1); % V = Vmax
% % % 
% % % 	    % Boundary condition for Smin and Smax
% % % 	    for v=1:NV-1
% % % 		    U(1,v) = 0;
% % % 		    U(NS,v) = max(0, Smax - K); % Galiotos uses U(NS-1,v) + ds;
% % % 	    end
% % % 	    % Boundary condition for Vmax
% % % 	    for s=1:NS
% % % 		    U(s,NV) = max(0, S(s) - K); % Galiotos uses U(s,NV-1);
% % % 	    end
% % % 	    % Update the temporary grid u(s,t) with the boundary conditions
% % % 	    u = U;
% % % 	    % Boundary condition for Vmin.
% % % 	    % Previous time step values are in the temporary grid u(s,t)
% % % 	    for s=2:NS-1
% % % 		    U(s,1) = u(s,1)*(1 - r*dt - kappa*theta*dt/dv)...
% % % 			       + dt*(0.5*(r-q)*s*u(s+1,1) - u(s-1,1)) ...
% % % 			       + kappa*theta*dt/dv*u(s,2);
% % % 	    end
% % % 
% % %         % Construct the tridiagonal matrix for each j
% % %         for j = 2:NV-1
% % %             A = zeros(NS-2, NS-2);
% % %             d = zeros(NS-2, 1);
% % % 
% % %             for i = 2:NS-1
% % % 
% % %                 alpha   = 0.25 * dt * (sigma^2 * V(j) * (S(i)^2) - kappa * (theta - V(j)) * S(i));
% % %                 beta    = 0.25 * dt * (sigma^2 * V(j) * (S(i)^2) + kappa * (theta - V(j)) * S(i));
% % %                 gamma   = 0.25 * dt * (r * S(i) - rho * sigma * V(j) * S(i) * j);
% % %                 delta   = 0.25 * dt * (r * S(i) + rho * sigma * V(j) * S(i) * j);
% % %                 epsilon =  0.5 * dt * (r * V(j) - rho * sigma * V(j) * j);
% % % 
% % %                 % Coefficients for the tridiagonal matrix
% % %                 a = -alpha - delta - epsilon;
% % %                 b = 1 + dt * r + 2 * epsilon;
% % %                 c = -beta - gamma - epsilon;
% % % 
% % %                 % Add a small value to diagonal to avoid singularity
% % %                 %b = b + 1e-10;
% % % 
% % %                 % Fill the tridiagonal matrix A
% % %                 if i > 2
% % %                     A(i-1, i-2) = a;
% % %                 end
% % %                 A(i-1, i-1) = b;  % This should be a scalar
% % %                 if i < NS-1
% % %                     A(i-1, i) = c;
% % %                 end
% % % 
% % %                 % Construct the right-hand side vector
% % %                 d(i-1) = alpha * U(i-1,j-1) + delta * U(i-1,j+1) + (1 - dt * r - 2 * epsilon) * U(i,j) + gamma * U(i+1,j-1) + beta * U(i+1,j+1);
% % %             end
% % % 
% % %             % Solve the linear system using a more robust method
% % %             u(2:NS-1, j) = A \ d;
% % %         end
% % % 
% % %         % Update the solution matrix
% % %         U = u;
% % % 
% % %         % % Enforce boundary conditions
% % %         % U(1, :) = 0; % S = 0
% % %         % U(NS, :) = Smax - K * exp(-r * (t * dt)); % S = Smax
% % %         % U(:, 1) = 0; % V = 0
% % %         % U(:, NV) = U(:, NV-1); % V = Vmax
% % %     end
% % % 
% % % end
% % % 

function U = HestonExplicitEuler(params,K,r,q,S,V,T,epsilon)

    % Finite differences for the Heston PDE for a European Call
    % Uses even grid sizes

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

    %single loop here producing X and Y
    for s=1:NS
	    for v=1:NV
		    U(s,v) = max(S(s) - K, 0);
	    end
    end

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
		    U(s,1) = u(s,1)*(1 - r*dt - kappa*theta*dt/dv)...
			       + dt*(0.5*(r-q)*s*u(s+1,1) - u(s-1,1)) ...
			       + kappa*theta*dt/dv*u(s,2);
	    end
	    % Update the temporary grid u(s,t) with the boundary conditions
	    u = U;

	    % Interior points of the grid (non boundary).
	    % Previous time step values are in the temporary grid u(s,t)

        % remember U indexation: U=U(s,v)
        % check page 317 in the book - this is explained
	    for s=2:NS-1
		    for v=2:NV-1
                %Here I will write the various components of the PDE
                derS=(u(s+1,v)-u(s-1,v))/(2*ds);
                derS2=(u(s+1,v)-2*u(s,v)+u(s-1,v))/(ds^2);
                derV=(u(s,v+1)-u(s,v-1))/(2*dv);
                derV2=(u(s,v+1)-2*u(s,v)+u(s,v-1))/(dv^2);
                derVS=(u(s+1,v+1)+u(s-1,v-1)-u(s-1,v+1)-u(s+1,v-1))/(4*dv*ds);

                sValue =S(s);
                vValue=V(v);

                first = (0.5*vValue*sValue^2*derS2);
                second = (0.5*sigma^2*vValue*derV2);
                third = (rho*sigma*vValue*sValue*derVS);
                fourth = ((r-q)*sValue*derS);
                fifth = (kappa*(theta-vValue)*derV)-r*u(s,v);

                U(s,v)=u(s,v)+dt*(first+second+third+fourth+fifth);
		    end
        end
    end
end




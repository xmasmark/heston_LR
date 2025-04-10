










function price = HestonVanillaClosedForm(S, K, T, r, kappa, theta, sigma, rho, v0, type)
    % Inputs:
    % S     = current stock price
    % K     = strike
    % T     = maturity
    % r     = risk-free rate
    % kappa = mean reversion speed
    % theta = long-run variance
    % sigma = volatility of variance
    % rho   = correlation between stock and volatility
    % v0    = current variance
    % type  = 'call' or 'put'
   
    % Integrals for P1 and P2
    P1 = 0.5 + (1/pi) * integral(@(phi) real(HestonIntegrandRL(phi, 1, S, K, T, r, kappa, theta, sigma, rho, v0)) ./ phi, 0, 100);
    P2 = 0.5 + (1/pi) * integral(@(phi) real(HestonIntegrandRL(phi, 2, S, K, T, r, kappa, theta, sigma, rho, v0)) ./ phi, 0, 100);
   
    % Call price from P1 and P2
    call = S * P1 - K * exp(-r*T) * P2;

    % Return call or put
    if strcmpi(type, 'call')
        price = call;
    elseif strcmpi(type, 'put')
        price = call - S + K * exp(-r*T);
    else
        error('type must be ''call'' or ''put''');
    end
end

function integrand = HestonIntegrandRL(phi, j, S, K, T, r, kappa, theta, sigma, rho, v0)
    % Characteristic function integrand
    i = complex(0,1);
   
    if j == 1
        u = 0.5;
        b = kappa - rho * sigma;
    else
        u = -0.5;
        b = kappa;
    end
   
    a = kappa * theta;
    x = log(S);
   
    d = sqrt((rho*sigma*i*phi - b).^2 - sigma^2*(2*u*i*phi - phi.^2));
    g = (b - rho*sigma*i*phi + d) ./ (b - rho*sigma*i*phi - d);
   
    C = (r*i*phi*T) + (a/sigma^2) * ((b - rho*sigma*i*phi + d)*T - 2*log((1 - g.*exp(d*T))./(1 - g)));
    D = ((b - rho*sigma*i*phi + d)/sigma^2) .* ((1 - exp(d*T))./(1 - g.*exp(d*T)));
   
    integrand = exp(C + D*v0 + i*phi*x) ./ (i*phi);
end


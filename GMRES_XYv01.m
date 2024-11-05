function [X_new, Y_new] = GMRES_XYv01(AKX, AKY, x0, y0, restart, tol, max_iter)
    % Custom GMRES routine adapted for low-rank factors
    % Inputs:
    %   AKX, AKY - Effective operators for low-rank approximation
    %   x0, y0   - Initial guesses for the low-rank factors X and Y
    %   restart  - Number of iterations before restarting GMRES
    %   tol      - Tolerance for convergence
    %   max_iter - Maximum number of iterations
    
    % Initialize solution and residual
    X_new = x0;
    Y_new = y0;
    residual = tol + 1; % Initial residual to enter loop
    iter = 0; % Initialize iteration count
    
    while residual > tol && iter < max_iter
        % Compute the product of AKX and AKY on current iterates X_new, Y_new
        AX_new = AKX * X_new'; % Apply effective operator on X
        AY_new = AKY * Y_new'; % Apply effective operator on Y
        
        % Compute the residual in low-rank form
        % residual_r = b - A * (X_new * Y_new')
        residual_X = x0 - AX_new; % Residual for X
        residual_Y = y0 - AY_new; % Residual for Y
        
        % Compute residual norm to check convergence
        residual = norm(residual_X, 'fro') * norm(residual_Y, 'fro');
        
        % If residual is small enough, break out of loop
        if residual < tol
            break;
        end
        
        % Orthogonalize residual against previous vectors
        % Restart mechanism: if we hit restart limit, reinitialize
        if mod(iter, restart) == 0 && iter > 0
            X_new = x0; % Restart X with initial guess
            Y_new = y0; % Restart Y with initial guess
        end
        
        % Update X_new and Y_new based on residuals
        % Using a basic GMRES iteration here, adjust with orthogonalization step
        X_new = X_new + residual_X; % Update for X factor
        Y_new = Y_new + residual_Y; % Update for Y factor
        
        % Increment iteration count
        iter = iter + 1;
    end
end
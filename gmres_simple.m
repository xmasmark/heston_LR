function [x, flag, relres, iter] = gmres_simple(A, b, tol, max_iter)
    % GMRES_SIMPLE - A basic implementation of the GMRES algorithm
    % Solves Ax = b using the GMRES method with a given tolerance and max iterations.
    %
    % Inputs:
    %   A        - Coefficient matrix (n x n)
    %   b        - Right-hand side vector (n x 1)
    %   tol      - Tolerance for convergence
    %   max_iter - Maximum number of iterations
    %
    % Outputs:
    %   x        - Solution vector
    %   flag     - Convergence flag (0 if converged, 1 if not)
    %   relres   - Relative residual norm
    %   iter     - Number of iterations performed

    n = length(b);
    x = zeros(n, 1); % Initial guess (zero vector)
    r = b - A * x;   % Initial residual
    beta = norm(r);
    
    if beta < tol
        flag = 0;
        relres = beta / norm(b);
        iter = 0;
        return;
    end

    V = zeros(n, max_iter + 1);
    H = zeros(max_iter + 1, max_iter);
    
    V(:,1) = r / beta;
    
    for k = 1:max_iter
        % Arnoldi process: Compute A*V_k and orthogonalize
        w = A * V(:,k);
        for j = 1:k
            H(j,k) = V(:,j)' * w;
            w = w - H(j,k) * V(:,j);
        end
        H(k+1,k) = norm(w);
        
        if H(k+1,k) ~= 0
            V(:,k+1) = w / H(k+1,k);
        end
        
        % Solve the least squares problem min || H*y - beta * e_1 ||
        e1 = zeros(k+1, 1);
        e1(1) = beta;
        y = H(1:k+1, 1:k) \ e1;
        
        % Compute the approximate solution
        x = V(:, 1:k) * y;
        
        % Compute residual norm
        relres = norm(A*x - b) / norm(b);
        
        if relres < tol
            flag = 0;
            iter = k;
            return;
        end
    end
    
    % If max iterations reached without convergence
    flag = 1;
    iter = max_iter;
end

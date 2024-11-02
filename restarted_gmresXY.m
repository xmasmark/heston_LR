function x = restarted_gmresXY(X, Y, b, x0, restart, tol, max_iter)
    % Restarted GMRES algorithm with low-rank approximation.
    % 
    % Parameters:
    % X, Y : matrices
    %     Low-rank approximation such that A ≈ X * Y'.
    % b : vector
    %     Right-hand side vector.
    % x0 : vector
    %     Initial guess for the solution.
    % restart : int
    %     Number of iterations before restart.
    % tol : float
    %     Tolerance for convergence.
    % max_iter : int
    %     Maximum number of iterations.
    %
    % Returns:
    % x : vector
    %     Solution vector.

    n = length(b);
    x = x0;
    r = b - X * (Y' * x);  % Adjust residual calculation
    beta = norm(r);
    
    for iter = 1:max_iter
        [Q, H] = arnoldi_processXY(X, Y, r, restart);  % Pass X, Y to arnoldi_process
        e1 = zeros(restart+1, 1);
        e1(1) = beta;
        
        % Solve the least squares problem
        [Q2, R] = qr(H,0);
        y = R \ (Q2' * e1);

        x = x + Q(:, 1:restart) * y;
        r = b - X * (Y' * x);  % Update residual calculation
        beta = norm(r);
        
        if beta < tol
            break;
        end
    end
end

function [Q, H] = arnoldi_processXY(X, Y, q, restart)
    % Arnoldi Process with low-rank approximation.
    % 
    % Parameters:
    % X, Y : matrices
    %     Low-rank approximation such that A ≈ X * Y'.
    % q : vector
    %     Initial vector for Krylov subspace.
    % restart : int
    %     Number of Arnoldi iterations.
    %
    % Returns:
    % Q : matrix
    %     Orthonormal basis for the Krylov subspace.
    % H : matrix
    %     Upper Hessenberg matrix.

    n = length(q);
    H = zeros(restart+1, restart);
    Q = zeros(n, restart+1);
    Q(:, 1) = q / norm(q);
    
    for k = 1:restart
        % Replace A * Q(:, k) with X * (Y' * Q(:, k))
        y = X * (Y' * Q(:, k));  
        for j = 1:k
            H(j, k) = Q(:, j)' * y;
            y = y - H(j, k) * Q(:, j);
        end
        H(k+1, k) = norm(y);
        if H(k+1, k) ~= 0 && k+1 <= restart
            Q(:, k+1) = y / H(k+1, k);
        end
    end
end

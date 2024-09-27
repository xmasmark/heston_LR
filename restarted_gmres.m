function x = restarted_gmres(A, b, x0, restart, tol, max_iter)
    % Restarted GMRES algorithm.
    % 
    % Parameters:
    % A : matrix
    %     Coefficient matrix.
    % b : vector
    %     Right-hand side vector.
    % x0 : vector
    %     Initial guess for the solution.
    % m : int
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
    r = b - A * x;
    beta = norm(r);
    
    for iter = 1:max_iter
        [Q, H] = arnoldi_process(A, r, restart);
        e1 = zeros(restart+1, 1);
        e1(1) = beta;
        
        % Solve the least squares problem
        %y = H \ e1; %this goes hysterical because of the sparsity

        [Q2, R] = qr(H,0);
        y = R \ (Q2'*e1);

        x = x + Q(:, 1:restart) * y;
        r = b - A * x;
        beta = norm(r);
        
        if beta < tol
            break;
        end
    end
end

function [Q, H] = arnoldi_process(A, q, restart)
    %as per definition, this implements the Arnoldi Process
    %meaning the matrix A, large and in this case sparse
    %is reduced to a Hessemberg matrix H with an orthonormal basis
    %Q for the Krylov subspace
    %please note that the Hessemberg matrix can be in upper or lower form
    %upper: a(i,j)~=0 for i,j with i>j+1;
    %lower: a(i,j)~=0 for i,j with j>i+1;
    %in both definitions, the +1 qualifies the defition that the Hessemberg
    %matrix is "almost" triangular.
    %REMEMBER that: AQ approx QH and this is an important proprerty that
    %makes this decomposition valuable in iterative algorithms.

    n = length(q);
    H = zeros(restart+1, restart);
    Q = zeros(n, restart+1);
    Q(:, 1) = q / norm(q);
    
    for k = 1:restart
        y = A * Q(:, k);
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

% function Hinv = invert_via_QR(A)
%     [Q,R] = qr(A);
%     % R_inv = inv(R); %inverse
%     % Q_T = Q'; %transposed
%     % Hinv = R_inv*Q_T;
% end
% 





















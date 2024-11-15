function [X_new, Y_new] = GMRES_XYv01(x,y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, K, Tmax, t, T, BX, BY, x0, y0, restart, tol, max_iter)


    %[BX,BY] = HestonMatVecBoundaries    (NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, K, Tmax, t, T)
    %function [xl,yl] = HestonMatVec(x,y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho)


   % Custom GMRES routine adapted for low-rank factors
    % Inputs:
    %   AKX, AKY - operators for low-rank approximation
    %   x0, y0   - Initial guesses for the low-rank factors X and Y
    %   restart  - Number of iterations before restarting GMRES
    %   tol      - Tolerance for convergence
    %   max_iter - Maximum number of iterations

    % Initialize solution and residual
    X_new = x0;
    Y_new = y0;
    % residual = tol + 1; % Initial residual to enter loop
    % iter = 0; % Initialize iteration count
    r = size(X_new, 2); % rank
    
    % Small threshold to avoid near-zero division issues
    % small_thresh = eps;  

    %r = b - A * x; x is replaced by X0 and Y0 
    %x is x0, y0
    %A*x is [xl,yl] = HestonMatVec(x0, y0, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho)
    %b is bx*by
    %r = BX*BY-xl*yl;
    %r is needed in low rank format
    %r = rx*ry;
    %r = BX*BY-xl*yl;
    %the factors need to be stacked-- [BX,+xl]*[BY,-yl] -- these are
    %stacked
    [xl,yl] = HestonMatVec(x,y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);
    r = BX*BY' - xl*yl'; %%b-A*X 

    first = [BX,xl];
    second = [BY,-yl];

    altR = first*second';
    beta = norm(r);

    for iter = 1:max_iter
        % [Q,H]= arnoldi_processXY(xl*yl', r, restart);
        %%[Q,H]= arnoldi_processXY(xl*yl', r(:), restart);
        [Q,H]= arnoldi_process_XY_I(xl*yl', r, restart);
        e1 = zeros(restart+1,1);
        e1(1)=beta;
        [Q2,R]=qr(H,0);
        y = R\(Q2'*e1);
        x = x+Q(:,1:restart)*y;
        r = b-A*x;
        beta = norm(r);

        if beta<tol
            break;
        end
    end

    % % n = length(b);
    % % x = x0;
    % % r = b - A * x;
    % % beta = norm(r);
    % % 
    % % for iter = 1:max_iter
    % %     [Q, H] = arnoldi_process(A, r, restart);
    % %     e1 = zeros(restart+1, 1);
    % %     e1(1) = beta;
    % % 
    % %     % Solve the least squares problem
    % %     %y = H \ e1; %this goes hysterical because of the sparsity
    % % 
    % %     [Q2, R] = qr(H,0);
    % %     y = R \ (Q2'*e1);
    % % 
    % %     x = x + Q(:, 1:restart) * y;
    % %     r = b - A * x;
    % %     beta = norm(r);
    % % 
    % %     if beta < tol
    % %         break;
    % %     end
    % % end



    
    

end


function [Q, H] = arnoldi_processXY(A, q, restart)
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

function [Q, H] = arnoldi_process_XY_I(A, q, restart)
    % arnoldi_process_flat: Arnoldi process with flattened A and q
    % A_flat: Flattened matrix A (A(:))
    % q_flat: Flattened starting vector q (q(:))
    % restart: Number of iterations
    % m: Number of rows in original matrix A
    % n: Number of columns in original matrix A


    % Ensure q is not zero
    if norm(q) == 0
        error('Starting vector q cannot be zero.');
    end

    % Initialize outputs
    H = zeros(restart + 1, restart); % Hessenberg matrix
    Q = zeros(length(q), restart + 1); % Krylov basis

    % Normalize q and set the first Krylov vector
    Q(:, 1) = q / norm(q);

    for k = 1:restart
        % Matrix-vector product: A * Q(:, k)
        y = A * Q(:, k);

        % Gram-Schmidt orthogonalization
        for j = 1:k
            H(j, k) = Q(:, j)' * y;
            y = y - H(j, k) * Q(:, j);
        end

        % Compute the next Hessenberg entry
        H(k + 1, k) = norm(y);

        % Break early if the norm is zero (linear dependence)
        if H(k + 1, k) == 0
            break;
        end

        % Normalize and add to Krylov basis
        if k + 1 <= restart
            Q(:, k + 1) = y / H(k + 1, k);
        end
    end
end

    % % n = length(b);
    % % x = x0;
    % % r = b - A * x;
    % % beta = norm(r);
    % % 
    % % for iter = 1:max_iter
    % %     [Q, H] = arnoldi_process(A, r, restart);
    % %     e1 = zeros(restart+1, 1);
    % %     e1(1) = beta;
    % % 
    % %     % Solve the least squares problem
    % %     %y = H \ e1; %this goes hysterical because of the sparsity
    % % 
    % %     [Q2, R] = qr(H,0);
    % %     y = R \ (Q2'*e1);
    % % 
    % %     x = x + Q(:, 1:restart) * y;
    % %     r = b - A * x;
    % %     beta = norm(r);
    % % 
    % %     if beta < tol
    % %         break;
    % %     end
    % % end






function [X_new, Y_new] = GMRES_LowRankV01(x,y, A, B, r, BX, BY, x0, y0, restart, tol, max_iter, dt)


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
    %rank = size(X_new, 2); % rank
    
    % Small threshold to avoid near-zero division issues
    % small_thresh = eps;  


    % n = length(b);
    % x = x0;
    % r = b - A * x;
    % beta = norm(r);
    % 


    %r = b - A * x; x is replaced by X0 and Y0 
    %A*x is [xl,yl] = HestonMatVec(x0, y0, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho)
    %b is bx*by
    %r = BX*BY-xl*yl;
    %r is needed in low rank format
    %r = rx*ry;
    %r = BX*BY-xl*yl;
    %the factors need to be stacked-- [BX,+xl]*[BY,-yl] -- these are
    %stacked
    %AX
    %[xl,yl] = HestonMatVec(x0,y0, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);
    [xl,yl] = LowRankMatVec(A,B,x0,y0);
    %residual = BX*BY' - xl*yl'; %%b-A*X 

    % residualX = [BX, -xl];
    % residualY = [BY, yl];


    residualX = [BX, -(1+dt*r/2)*x0, (dt/2)*xl];
    residualY = [BY,             y0,         yl];

    beta = norm_lr(residualX,residualY);

    for iter = 1:max_iter
        
        %[Qx, Qy, H]= arnoldi_process_low_rank_steroids(residualX, residualY, A, B, restart, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, tol, dt);
        [Qx, Qy, H]= arnoldi_process_low_rank_steroids(residualX, residualY, A, B, restart, r, tol, dt);
        e1 = zeros(restart+1,1);
        e1(1)=beta;
        [Q2,R]=qr(H,0);
        yA = R\(Q2'*e1);
        %next becomes a low rank process
        %loop 1 - restart, at each step I add the vector of the cellarray
        %to the low rank factor via Qx and Qy
        for i = 1:restart
           %x = x+Q(:,r)*yA(r);
           x = [x, Qx{i}*yA(i)];
           y = [y, Qy{i}];
        end
        %how about compressing x and y now before they go further into the
        %program?
        %x = x+Q(:,1:restart)*y;
        %r = b-A*x;
        %[Ax,Ay] = HestonMatVec(x,y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);
        [Ax,Ay] = LowRankMatVec(A,B,x,y);
        % residualX = [BX, -Ax];
        % residualY = [BY, Ay];

        residualX = [BX, -(1+dt*r/2)*x, (dt/2)*Ax];
        residualY = [BY,             y,        Ay];
        

        [residualX,residualY] = CompressData(residualX,residualY,tol);

        %beta recalculation
        beta = norm_lr(residualX,residualY);

        X_new = x;
        Y_new = y;

        if beta<tol
            break;
        end
    end
end


%[Q,H]= arnoldi_process_XY_I(xl, yl, residualX, residualY, restart);
%function [Qx, Qy, H] = arnoldi_process_low_rank_steroids(residualX, residualY, A, B, restart, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, tol, dt)
function [Qx, Qy, H] = arnoldi_process_low_rank_steroids(residualX, residualY, A, B, restart, r, tol, dt)
    % arnoldi_process_flat: Arnoldi process with flattened A and q
    % A_flat: Flattened matrix A (A(:))
    % q_flat: Flattened starting vector q (q(:))
    % restart: Number of iterations
    % m: Number of rows in original matrix A
    % n: Number of columns in original matrix A


    % % Ensure q is not zero
    % if norm(q) == 0
    %     error('Starting vector q cannot be zero.');
    % end

    %replacement for Q:
    Qx = cell(restart+1,1);
    Qy = cell(restart+1,1);

    % Initialize outputs
    H = zeros(restart + 1, restart); % Hessenberg matrix
    %Q = zeros(length(q), restart + 1); % Krylov basis
    %Q = zeros(length(q), restart + 1); % Krylov basis


    % Normalize q and set the first Krylov vector
    % is meant to be a vector and it's a matrix instead
    % so that raises an error in the following statement:
    %Q(:, 1) = q / norm(q);

    %The Norm needs to be computed as discussed in low rank format
    normXY = norm_lr(residualX,residualY);
    % Qx{1}=residualX/norm(residualX);
    % Qy{1}=residualY/norm(residualY);

    Qx{1}=residualX/(normXY^0.5);
    Qy{1}=residualY/(normXY^0.5);

    for k = 1:restart
        % Matrix-vector product: A * Q(:, k)
        %y = A * Q(:, k);

        %this is the equivalent of A*Q
        %[Yx,Yy] = HestonMatVec(Qx{k},Qy{k}, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);
        [Yx,Yy] = LowRankMatVec(A,B,Qx{k},Qy{k});

        %change of sign because not residuals
        Yx = [(1+dt*r/2)*Qx{k}, -(dt/2)*Yx];
        Yy = [Qy{k}, Yy];

        % Gram-Schmidt orthogonalization
        for j = 1:k
            [Qxj,Qyj] = CompressData(Qx{j}, Qy{j}, tol);
            H(j, k) = dot_lr(Qxj,Qyj,Yx,Yy);

            %y = y - H(j, k) * Q(:, j);

            Yx = [Yx, H(j,k)*Qx{j}];
            Yy = [Yy, -Qy{j}];

            [Yx,Yy] = CompressData(Yx,Yy,tol);
        end

        % Compute the next Hessenberg entry
        % Mind the norm -- already discussed and documented
        %H(k + 1, k) = norm(y);
        
        normv = norm_lr(Yx,Yy);
        H(k + 1, k) = normv;

        % Break early if the norm is zero (linear dependence)
        if H(k + 1, k) == 0
            break;
        end

        % Normalize and add to Krylov basis
        if k + 1 <= restart
            %Q(:, k + 1) = y / H(k + 1, k);
            Qx{k+1}=Yx/H(k + 1, k);
            Qy{k+1}=Yy; %I have to divide only once
        end
    end
end



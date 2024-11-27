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
    [xl,yl] = HestonMatVec(x0,y0, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);
    %residual = BX*BY' - xl*yl'; %%b-A*X 

    residualX = [BX, -xl];
    residualY = [BY, yl];

    beta = norm_lr(residualX,residualY);

    % residual_compressed = sum(residual, 2);
    % 
    % beta = norm(residual);

    %the issue is that residual is not a vector but a NSxNV matrix!
    %so the arnoldi_Process_XY_I cannot work:

    for iter = 1:max_iter
        % [Q,H]= arnoldi_processXY(xl*yl', r, restart);
        %%[Q,H]= arnoldi_processXY(xl*yl', r(:), restart);
        %[Q,H]= arnoldi_process_XY_I(xl*yl', residual, restart);
        % [Q,H]= arnoldi_process_XY_I(xl, yl, residualX, residualY, restart);
        [Qx, Qy, H]= arnoldi_process_low_rank(residualX, residualY, restart, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, tol);
        e1 = zeros(restart+1,1);
        e1(1)=beta;
        [Q2,R]=qr(H,0);
        yA = R\(Q2'*e1);
        %next becomes a low rank process
        %loop 1 - restart, at each step I add the vector of the cellarray
        %to the low rank factor via Qx and Qy
        for r = 1:restart
           %x = x+Q(:,r)*yA(r);
           x = [x, Qx{r}*yA(r)];
           y = [y, Qy{r}];
        end
        %how about compressing x and y now before they go further into the
        %program?
        %x = x+Q(:,1:restart)*y;
        %r = b-A*x;
        [Ax,Ay] = HestonMatVec(x,y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);
        residualX = [BX, -Ax];
        residualY = [BY, Ay];

        [residualX,residualY] = CompressData(residualX,residualY,tol);

        %new norm that I will implement taking as parameters residualX,
        %residualY
        % beta = norm(r);
        beta = norm_lr(residualX,residualY);

        if beta<tol
            break;
        end
        X_new = x;
        Y_new = y;
    end
end


% % % function [Q, H] = arnoldi_processXY(A, q, restart)
% % %     %as per definition, this implements the Arnoldi Process
% % %     %meaning the matrix A, large and in this case sparse
% % %     %is reduced to a Hessemberg matrix H with an orthonormal basis
% % %     %Q for the Krylov subspace
% % %     %please note that the Hessemberg matrix can be in upper or lower form
% % %     %upper: a(i,j)~=0 for i,j with i>j+1;
% % %     %lower: a(i,j)~=0 for i,j with j>i+1;
% % %     %in both definitions, the +1 qualifies the defition that the Hessemberg
% % %     %matrix is "almost" triangular.
% % %     %REMEMBER that: AQ approx QH and this is an important proprerty that
% % %     %makes this decomposition valuable in iterative algorithms.
% % % 
% % %     n = length(q);
% % %     H = zeros(restart+1, restart);
% % %     Q = zeros(n, restart+1);
% % %     Q(:, 1) = q / norm(q);
% % % 
% % %     for k = 1:restart
% % %         y = A * Q(:, k);
% % %         for j = 1:k
% % %             H(j, k) = Q(:, j)' * y;
% % %             y = y - H(j, k) * Q(:, j);
% % %         end
% % %         H(k+1, k) = norm(y);
% % %         if H(k+1, k) ~= 0 && k+1 <= restart
% % %             Q(:, k+1) = y / H(k+1, k);
% % %         end
% % %     end
% % % end

%[Q,H]= arnoldi_process_XY_I(xl, yl, residualX, residualY, restart);
function [Qx, Qy, H] = arnoldi_process_low_rank(residualX, residualY, restart, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, tol)
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

    Qx{1}=residualX/normXY;
    Qy{1}=residualY/normXY;

    for k = 1:restart
        % Matrix-vector product: A * Q(:, k)
        %y = A * Q(:, k);

        [Yx,Yy] = HestonMatVec(Qx{k},Qy{k}, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);

        % Gram-Schmidt orthogonalization
        for j = 1:k
            %this is calculated with the custom dot product being careful
            %in the factors rearrangement because if I don't do that, I
            %lose the low rank efficiency
            % H(j, k) = Q(:, j)' * y;
            % H(j, k) = dot_lr(Qx{},Qy{});

            %two parameters... not four... but they need to be compressed
            %or the thing will blow up
            [Qxj,Qyj] = CompressData(Qx{j},Qy{j},tol);
            [Qxk,Qyk] = CompressData(Qx{k},Qy{k},tol);
            % H(j, k) = dot_lr(Qx{j},Qy{j},Qx{k},Qy{k});
            H(j, k) = dot_lr(Qxj,Qyj,Qxk,Qyk);

            %y = y - H(j, k) * Q(:, j);

            %y will be a low rank vector calculated 
            %the result will be two new Yx and Yy calculated from the above
            %results
            Yx = [Yx, H(j,k)*Qx{k}];
            Yy = [Yy, -Qy{k}];

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



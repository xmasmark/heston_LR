function [X_new, Y_new] = GMRES_XY_COMP(A, b, x0v, x, y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, K, Tmax, t, T, BX, BY, x0, y0, restart, tol, max_iter, dt)

    testb = reshape(b,[NS,NV]) - BX*BY';

    fnb = norm(testb,'fro');

    testInitial = reshape(x0v,[NS,NV])-x0*y0';
    fnIni = norm(testInitial,'fro');

    % Initialize solution and residual
    X_new = x0;
    Y_new = y0;

    [xl,yl] = HestonMatVec(x0, y0, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);

    %this is the matrix multiplied by A however, it's not what I need for
    %the linear system..
    %lhs_matrix = ((1+(dt*r/2))*eye(NS*NV) - (dt/2)*AK);  
    %lhs matrix * x0v = ((1+(dt*r/2))*eye(NS*NV)*x0v - (dt/2)*AK*x0v);
    %lhs matrix * x0v = ((1+(dt*r/2))*x0v - (dt/2)*xlv);
    %lhs matrix * x0v = (1+dt*r/2)*x0v - dt/2*xlv;
    % residualX = [BX, -xl*dt];
    % residualY = [BY, yl];

    residualX = [BX, -(1+dt*r/2)*x0, (dt/2)*xl];
    residualY = [BY,             y0,         yl];
    
    beta = norm_lr(residualX,residualY);

    for iter = 1:max_iter
        
        [Qx, Qy, Hlr]= arnoldi_process_low_rank(residualX, residualY, restart, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, tol);
        e1lr = zeros(restart+1,1);
        e1lr(1)=beta;

        % [Q2lr,Rlr]=qr(Hlr,0);
        % yA = Rlr\(Q2lr'*e1lr);

        yA=Hlr\e1lr;

        %next becomes a low rank process
        %loop 1 - restart, at each step I add the vector of the cellarray
        %to the low rank factor via Qx and Qy
        for i = 1:restart
           %x = x+Q(:,r)*yA(r);
           x = [x, Qx{i}*yA(i)];
           y = [y, Qy{i}];
        end

        %x = x+Q(:,1:restart)*y;
        %r = b-A*x;

        [Ax,Ay] = HestonMatVec(x, y, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);
        % residualX = [BX, -Ax];
        % residualY = [BY,  Ay];

        residualX = [BX, -(1+dt*r/2)*x, (dt/2)*Ax];
        residualY = [BY,             y,        Ay];
        

        [residualX,residualY] = CompressData(residualX, residualY, tol);

        %new norm that I will implement taking as parameters residualX,
        %residualY
        % beta = norm(r);
        beta = norm_lr(residualX, residualY);

        X_new = x;
        Y_new = y;
        
        if beta<tol
            break;
        end
        %result of low rank
    end

    %n = length(b);
    xV = x0v;
    residual = b - A * xV;
    betaV = norm(residual);

    for iter = 1:max_iter
        [Q, H] = arnoldi_process(A, residual, restart);
        e1 = zeros(restart+1, 1);
        e1(1) = betaV;
        
        % Solve the least squares problem
        %y = H \ e1; %this goes hysterical because of the sparsity

        [Q2, R] = qr(H,0);
        yv = R \ (Q2'*e1);

        xV = xV + Q(:, 1:restart) * yv;
        residual = b - A * xV;
        betaV = norm(residual);
        
        U = reshape(xV, [NS, NV]);

        % if betaV < tol || beta < tol
        %     break;
        % end
        if beta < tol
            break;
        end
    end

    diff = U - X_new*Y_new';
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


%[Q,H]= arnoldi_process_XY_I(xl, yl, residualX, residualY, restart);
function [Qx, Qy, H] = arnoldi_process_low_rank(residualX, residualY, restart, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho, tol)

    %replacement for Q:
    Qx = cell(restart+1,1);
    Qy = cell(restart+1,1);

    % Initialize outputs
    H = zeros(restart + 1, restart); % Hessenberg matrix

    %The Norm needs to be computed as discussed in low rank format
    normXY = norm_lr(residualX,residualY);
    % Qx{1}=residualX/norm(residualX);
    % Qy{1}=residualY/norm(residualY);

    nr = sqrt(normXY);

    Qx{1}=residualX/(nr);
    Qy{1}=residualY/(nr);

    for k = 1:restart
        % Matrix-vector product: A * Q(:, k)
        %y = A * Q(:, k);

        [Yx,Yy] = HestonMatVec(Qx{k},Qy{k}, NS, NV, ds, dv, S, V, r, q, kappa, theta, lambda, sigma, rho);

        % 
        % residualX = [BX, -(1+dt*r/2)*x, (dt/2)*Ax];
        % residualY = [BY,             y,        Ay];

        %change of sign because not residuals
        Yx = [(1+dt*r/2)*Qx{k}, -(dt/2)*Yx];
        Yy = [Qy{k}, Yy];

        % Gram-Schmidt orthogonalization
        for j = 1:k
            [Qxj,Qyj] = CompressData(Qx{j},Qy{j},tol);
            % [Qxk,Qyk] = CompressData(Qx{k},Qy{k},tol);
            %H(j, k) = dot_lr(Qxj,Qyj,Qxk,Qyk);
            H(j, k) = dot_lr(Qxj,Qyj,Yx,Yy);

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



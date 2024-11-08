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
    r = size(X_new, 2); % rank
    
    % Loop until convergence or maximum iterations
    while residual > tol && iter < max_iter
        % Initialize Krylov subspace for each rank component
        V_X = zeros(size(X_new, 1), restart, r);
        V_Y = zeros(size(Y_new, 1), restart, r);
        H = zeros(restart + 1, restart, r);

        % Calculate initial residual and its norm for each rank component
        for k = 1:r

            % % % ISSUE HERE:
            % % % AKX is NSxr
            % % % X_new(:, k) is NSx1
            % % % 
            % % % AKY is NVxr
            % % % Y_new(:, k) is NVx1

            % AX_new = AKX * X_new(:, k); % Apply effective operator on each column of X
            % AY_new = AKY * Y_new(:, k); % Apply effective operator on each column of Y
            
            AX_new = AKX(:, k) .* X_new(:, k);
            AY_new = AKY(:, k) .* Y_new(:, k);


            residual_X = x0(:, k) - AX_new;
            residual_Y = y0(:, k) - AY_new;
            residual_k = norm(residual_X, 'fro') * norm(residual_Y, 'fro');

            % Normalize initial residual
            V_X(:, 1, k) = residual_X / norm(residual_X, 'fro');
            V_Y(:, 1, k) = residual_Y / norm(residual_Y, 'fro');
        end
        
        % Check total residual across rank components
        residual = sum(arrayfun(@(k) norm(V_X(:, 1, k), 'fro') * norm(V_Y(:, 1, k), 'fro'), 1:r));
        [num_rows, ~] = size(H(:, :, k));  % Dynamically get the number of rows of H
        residual_vector = norm(V_X(:, 1, k), 'fro') * norm(V_Y(:, 1, k), 'fro') * ones(num_rows, 1);

        % Arnoldi process to build Krylov subspace
        for j = 1:restart
            for k = 1:r
                % Apply the operators AKX and AKY on each rank component
                %PROBLEM HERE:
                %AKX is (NS,r) and V_X(:,j,k) is (NS,1) so the product
                %between the two cannot work
                w_X = AKX .* V_X(:, j, k);
                w_Y = AKY .* V_Y(:, j, k);



                % w_X = AKX' * V_X(:, j, k);
                % w_Y = AKY' * V_Y(:, j, k);


                % Orthogonalize w against previous V's
                for i = 1:j
                    %H(i, j, k) = V_X(:, i, k)' * w_X + V_Y(:, i, k)' * w_Y;
                    %H(i, j, k) = sum(V_X(:, i, k) .* w_X) + sum(V_Y(:, i, k) .* w_Y);
                    % H(i, j, k) = dot(V_X(:, i, k), w_X) + dot(V_Y(:, i, k), w_Y);
                    H(i, j, k) = sum(sum(V_X(:, i, k) .* w_X) + sum(V_Y(:, i, k) .* w_Y));

                    w_X = w_X - H(i, j, k) * V_X(:, i, k);
                    w_Y = w_Y - H(i, j, k) * V_Y(:, i, k);
                end
                H(j + 1, j, k) = norm(w_X, 'fro') * norm(w_Y, 'fro');

                % Normalize and add new basis vectors if H(j+1, j) is not zero
                if H(j + 1, j, k) > tol
                    V_X(:, j + 1, k) = w_X / norm(w_X, 'fro');
                    V_Y(:, j + 1, k) = w_Y / norm(w_Y, 'fro');
                else
                    break;
                end
            end
        end

        % Solve least-squares problem for updating X and Y for each rank component
        for k = 1:r
            y_k = H(:, :, k) \ residual_vector; % Solve H*y = residual for each rank component
            X_update = V_X(:, 1:j, k) * y_k(1:j); % Update for X
            Y_update = V_Y(:, 1:j, k) * y_k(1:j); % Update for Y

            % Update X_new and Y_new with computed update
            X_new(:, k) = X_new(:, k) + X_update;
            Y_new(:, k) = Y_new(:, k) + Y_update;
        end

        % Increment iteration count
        iter = iter + 1;
        
        % Restart if needed
        if mod(iter, restart) == 0
            X_new = x0;
            Y_new = y0;
        end
    end
end



% % % function [X_new, Y_new] = GMRES_XYv01(AKX, AKY, x0, y0, restart, tol, max_iter)
% % %     % Custom GMRES routine adapted for low-rank factors
% % %     % Inputs:
% % %     %   AKX, AKY - Effective operators for low-rank approximation
% % %     %   x0, y0   - Initial guesses for the low-rank factors X and Y
% % %     %   restart  - Number of iterations before restarting GMRES
% % %     %   tol      - Tolerance for convergence
% % %     %   max_iter - Maximum number of iterations
% % % 
% % %     % Initialize solution and residual
% % %     X_new = x0;
% % %     Y_new = y0;
% % %     residual = tol + 1; % Initial residual to enter loop
% % %     iter = 0; % Initialize iteration count
% % % 
% % %     % Loop until convergence or maximum iterations
% % %     while residual > tol && iter < max_iter
% % %         % Initialize Krylov subspace
% % %         % V_X = zeros(size(X_new, 1), restart);
% % %         % V_Y = zeros(size(Y_new, 1), restart);
% % %         V_X = zeros(size(X_new, 1), restart);
% % %         V_Y = zeros(size(Y_new, 1), restart);
% % %         H = zeros(restart + 1, restart);
% % % 
% % %         % Calculate initial residual and its norm
% % %         AX_new = AKX * X_new'; % Apply effective operator on X
% % %         AY_new = AKY * Y_new'; % Apply effective operator on Y
% % % 
% % %         % ISSUE: 
% % %         % AX_new is NSxNS
% % %         % AY_new is NVxNV
% % %         % x0 is NSxr
% % %         % y0 is NVxr
% % % 
% % %         residual_X = x0 - AX_new;
% % %         residual_Y = y0 - AY_new;
% % %         residual = norm(residual_X, 'fro') * norm(residual_Y, 'fro');
% % % 
% % %         % Normalize initial residual
% % %         V_X(:, 1) = residual_X / norm(residual_X, 'fro');
% % %         V_Y(:, 1) = residual_Y / norm(residual_Y, 'fro');
% % % 
% % %         % Arnoldi process to build Krylov subspace
% % %         for j = 1:restart
% % %             % Apply the operators AKX and AKY
% % %             w_X = AKX * V_X(:, j);
% % %             w_Y = AKY * V_Y(:, j);
% % % 
% % %             % Orthogonalize w against previous V's
% % %             for i = 1:j
% % %                 H(i, j) = V_X(:, i)' * w_X + V_Y(:, i)' * w_Y;
% % %                 w_X = w_X - H(i, j) * V_X(:, i);
% % %                 w_Y = w_Y - H(i, j) * V_Y(:, i);
% % %             end
% % %             H(j + 1, j) = norm(w_X, 'fro') * norm(w_Y, 'fro');
% % % 
% % %             % Normalize and add new basis vectors if H(j+1, j) is not zero
% % %             if H(j + 1, j) > tol
% % %                 V_X(:, j + 1) = w_X / norm(w_X, 'fro');
% % %                 V_Y(:, j + 1) = w_Y / norm(w_Y, 'fro');
% % %             else
% % %                 break;
% % %             end
% % %         end
% % % 
% % %         % Solve least-squares problem for updating X and Y
% % %         y = H \ residual; % Solve H*y = residual
% % %         X_update = V_X(:, 1:j) * y(1:j); % Update for X
% % %         Y_update = V_Y(:, 1:j) * y(1:j); % Update for Y
% % % 
% % %         % Update X_new and Y_new with computed update
% % %         X_new = X_new + X_update;
% % %         Y_new = Y_new + Y_update;
% % % 
% % %         % Recalculate residual for convergence check
% % %         residual_X = x0 - (AKX * X_new');
% % %         residual_Y = y0 - (AKY * Y_new');
% % %         residual = norm(residual_X, 'fro') * norm(residual_Y, 'fro');
% % % 
% % %         % Increment iteration count
% % %         iter = iter + 1;
% % % 
% % %         % Restart if needed
% % %         if mod(iter, restart) == 0
% % %             X_new = x0;
% % %             Y_new = y0;
% % %         end
% % %     end
% % % end

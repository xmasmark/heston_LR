% function [X, Y] = CompressData(XC, YC, epsilon)
%     % Economical QR decompositions
%     [Qx, Rx] = qr(XC, "econ");
%     [Qy, Ry] = qr(YC, "econ");
% 
%     % Intermediate product
%     Ri = Rx * Ry';
% 
%     % Compute dominant SVD and rank using faster method
%     [r, u, s, v] = SVDRankPlus(Ri, epsilon);
% 
%     % Extract dominant components
%     du = u(:, 1:r);
%     ds = s(1:r, 1:r);
%     dv = v(:, 1:r);
% 
%     % Precompute sqrt of ds once
%     Ssqrt = sqrt(ds);
% 
%     % Build compressed matrices
%     X = Qx * du * Ssqrt;
%     Y = Qy * dv * Ssqrt;
% end
% 
% function [R, uc, sc, vs] = SVDRankPlus(U, epsilon)
%     % Compute economy SVD only as needed
%     [uc, sc, vs] = svd(U, "econ");
% 
%     % Vector of squared singular values
%     svals = diag(sc).^2;
%     total = sum(svals);
% 
%     % Compute rank by cumulative sum from the smallest singular values
%     cutoffs = cumsum(flip(svals));
%     index = find(cutoffs > epsilon^2 * total, 1, 'first');
% 
%     % Reconstruct rank from flipped index
%     R = length(svals) - index + 2;
% 
%     % Safety clamp
%     R = min(R, length(svals));
% end

%old code

function [X, Y]=CompressData(XC,YC, epsilon)

    %[Qx,Rx]=qr(XC,0,"econ");
    %[Qy,Ry]=qr(YC,0,"econ");

    [Qx,Rx]=qr(XC,"econ");
    [Qy,Ry]=qr(YC,"econ");

    Ri = Rx*Ry';

    [rank, u, s, v]=SVDRankPlus2(Ri,epsilon);

    %the prefix d stands for "dominant"
    du = u(:,1:rank);
    ds = s(1:rank,1:rank);
    %dv = v(1:rank,:);   
    dv = v(:,1:rank);   

    %creation of compressed matrices
    X=Qx*du*sqrt(ds);
    Y=Qy*dv*sqrt(ds);

end


function [R, uc, sc, vs] = SVDRankPlus(U, epsilon)

    %builtin_rank = rank(U,epsilon);

    [uc,sc,vs]=svd(U, "econ");

    sigSum = 0;

    max=size(sc,2);
    if(size(sc,1)<size(sc,2))
        max = size(sc,1);
    end

    for i = 1: max
        sigSum = sigSum + power(sc(i,i),2);
    end

    threshold = power(epsilon,2)*sigSum;

    rank_home_made = 0;

    errorSquare = 0;

    for i=size(sc,2):-1:1
        errorSquare = errorSquare + power(sc(i,i),2);
        if errorSquare > threshold
            rank_home_made = i + 2;
            break;
        end
    end

    if(rank_home_made>max)
        rank_home_made = max;
    end

    R=rank_home_made;

end

%function [R, uc, sc, vs, info] = SVDRankPlus2(U, epsilon)
function [R, uc, sc, vs] = SVDRankPlus2(U, epsilon)
    [uc, sc, vs] = svd(U, "econ");
    svals = diag(sc).^2;
    sigSum = sum(svals);
    threshold = epsilon^2 * sigSum;

    % Accumulate tail until above threshold
    errorSquare = 0;
    R = 0;
    for i = length(svals):-1:1
        errorSquare = errorSquare + svals(i);
        if errorSquare > threshold
            R = i + 2;
            break;
        end
    end
    R = min(R, min(size(sc)));  % clamp
    info.svals = svals;
end
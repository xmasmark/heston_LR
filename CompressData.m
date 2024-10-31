function [X, Y]=CompressData(XC,YC, epsilon)

    %[Qx,Rx]=qr(XC,0,"econ");
    %[Qy,Ry]=qr(YC,0,"econ");

    [Qx,Rx]=qr(XC,"econ");
    [Qy,Ry]=qr(YC,"econ");

    Ri = Rx*Ry';
    
    [rank, u, s, v]=SVDRankPlus(Ri,epsilon);

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


function [xl,yl] = LowRankMatVec(A,B,x,y)

    szA = size(A);
    szB = size(B);

    xl = [];
    yl = [];

    for n = 1:szA(3)
        xl=[xl,A(:,:,n)*x];
    end

    for n = 1:szB(3)
        yl=[yl,B(:,:,n)*y];
    end

    % 
    % X1 = A1KX*x; %first slice of A A(:,:,1)
    % Y1 = A1KY*y; %first slice of B B(:,:,1)
    % X2 = A2KX*x;
    % Y2 = A2KY*y;
    % X3 = A3KX*x;
    % Y3 = A3KY*y;
    % X4 = A4KX*x;
    % Y4 = A4KY*y;
    % X5 = A5KX*x;
    % Y5 = A5KY*y;
    % 
    % xl = [X1,X2,X3,X4,X5];
    % yl = [Y1,Y2,Y3,Y4,Y5];
end
function MU = HestonExplicitLowRankRTildaRC2(params,K,r,q,S,V,T,epsilon)

% Finite differences for the Heston PDE for a European Call
% Uses even grid sizes
% In 'T Hout and Foulon "ADI Finite Difference Schemes for Option Pricing
% in the Heston Modelo with Correlation" Int J of Num Analysis and Modeling, 2010.
% Thesis by Sensi Li and paper by Vassilis Galiotos
% INPUTS
%    params = 6x1 vector of Heston parameters
%    K = Strike price
%    r = risk free rate
%    q = Dividend yield
%    S = vector for stock price grid
%    V = vector for volatility grid
%    T = vector for maturity grid
% OUTPUT
%    U = U(S,v) 2-D array of size (nS+1)x(nV+1) for the call price

% Heston parameters
kappa = params(1);
theta = params(2);
sigma = params(3);
v0    = params(4);
rho   = params(5);
lambda = params(6);
%epsilon = 0.001;

% Length of stock price, volatility, and maturity
NS = length(S);
NV = length(V);
NT = length(T);
Smin = S(1);  Smax = S(NS);
Vmin = V(1);  Vmax = V(NV);
Tmin = T(1);  Tmax = T(NT);

% Increment for Stock Price, Volatility, and Maturity
ds = (Smax-Smin)/(NS-1);
dv = (Vmax-Vmin)/(NV-1);
dt = (Tmax-Tmin)/(NT-1);

% Initialize the 2-D grid with zeros
U = zeros(NS,NV);

% Temporary grid for previous time steps
u = zeros(NS,NV);

% Solve the PDE
% Round each value of U(S,v,t) at each step
% Boundary condition for t = Maturity

%single loop here producing X and Y
for s=1:NS
	for v=1:NV
		U(s,v) = max(S(s) - K, 0);
	end
end

X=zeros(NS,1);
for s=1:NS
    X(s)=max(S(s)-K,0);
end

Y=zeros(NV,1);
for v=1:NV
	Y(v)=1;
end

%vectorised U = X*Y' where apostrophe is 
%the sign for transposition (column to row).

%testMatrix = U-X*Y';

    for t=1:NT-1
	    % Boundary condition for Smin and Smax
	    for v=1:NV-1
		    U(1,v) = 0;
		    %U(NS,v) = max(0, Smax - K); % Galiotos uses U(NS-1,v) + ds;
            U(NS,v)= U(NS-1,v) + ds;
	    end
	    % Boundary condition for Vmax
	    for s=1:NS
		    %U(s,NV) = max(0, S(s) - K); % Galiotos uses U(s,NV-1);
            %U(s,NV) = U(s,NV-1);
            U(s,NV) = S(s);
	    end
	    % Update the temporary grid u(s,t) with the boundary conditions
	    u = U;
	    % Boundary condition for Vmin.
	    % Previous time step values are in the temporary grid u(s,t)
	    for s=2:NS-1
		    %DerV = (u(s,2) - u(s,1)) / dv;                	% PDE Points on the middle of the grid (non boundary)
		    %DerS = (u(s+1,1) - u(s-1,1))/2/ds;              % Central difference for dU/dS
		    U(s,1) = u(s,1)*(1 - r*dt - kappa*theta*dt/dv)...
			       + dt*(0.5*(r-q)*s*u(s+1,1) - u(s-1,1)) ...
			       + kappa*theta*dt/dv*u(s,2);
	    end
	    % Update the temporary grid u(s,t) with the boundary conditions
	    u = U;
    
	    % Interior points of the grid (non boundary).
	    % Previous time step values are in the temporary grid u(s,t)
    
        % remember U indexation: U=U(s,v)
        % check page 317 in the book - this is explained
        copyOfU=U;
	    for s=2:NS-1
		    for v=2:NV-1
                %Here I will write the various components of the PDE
                derS=(u(s+1,v)-u(s-1,v))/(2*ds);
                derS2=(u(s+1,v)-2*u(s,v)+u(s-1,v))/(ds^2);
                derV=(u(s,v+1)-u(s,v-1))/(2*dv);
                derV2=(u(s,v+1)-2*u(s,v)+u(s,v-1))/(dv^2);
                derVS=(u(s+1,v+1)+u(s-1,v-1)-u(s-1,v+1)-u(s+1,v-1))/(4*dv*ds);
    
                sValue =S(s);
                vValue=V(v);
    
                first = (0.5*vValue*sValue^2*derS2);
                second = (0.5*sigma^2*vValue*derV2);
                third = (sigma*vValue*sValue*derVS);
                fourth = ((r-q)*sValue*derS);
                fifth = (kappa*(theta-vValue)*derV);
                sixth = -r*u(s,v);
    
                U(s,v)=u(s,v)+dt*(first+second+third+fourth+fifth+sixth);
                %U(s,v)=u(s,v)+dt*(first+third+fourth+fifth+sixth);
		    end
        end

        xo=X;
        yo=Y;

        %Boundary Conditions start

        xones = ones(NS,1);
        xzeros = zeros(NS,1);
        xzeros2 = zeros(NS,1);
        yzeros = zeros(NV,1);
        yzeros2 = zeros(NV,1);
        yzeros(NV)=1;
        yzeros2(1)=1;

        x = [X xzeros2 xzeros];
        y = [Y yzeros yzeros2];

        %first step is to 

        sz = size(x);
        records = sz(1);
        fields = sz(2);


        %U(s,NV) = S(s);

		    % U(s,1) = u(s,1)*(1 - r*dt - kappa*theta*dt/dv)+dt*(0.5*(r-q)*s*u(s+1,1) - u(s-1,1))+ kappa*theta*dt/dv*u(s,2);
            % U(s,1) = us1*(1 - r*dt - kappa*theta*dt/dv)+dt*(0.5*(r-q)*s*u(s+1,1) - u(s-1,1))+ kappa*theta*dt/dv*u(s,2);

        %U(s,1) = ...;
        %boundary condition for the first column of the matrix U
        %once this part of the code will reach maturity
        %both this loop and the previous one will be merged into one

        yc=Y(1,:);
        yc2=Y(2,:);
        ynv = Y(NV,:);

        for j = 2:(records-1)
            xc=X(j,:);
            xcp1=X(j+1,:);
            xcm1=X(j-1,:);

            us1 = xc*yc';
            usp1= xcp1*yc';
            us2 = xc*yc2';
            usm1= xcm1*yc';
            
            %new value, hence the n prefix:
            nus1=us1*(1-r*dt-kappa*theta*dt/dv)+dt*(0.5*(r-q)*j*usp1-usm1)+kappa*theta*dt/dv*us2;

            %acc=0;
            %in acc I need to store the cross 
            %product between X and the first column of Y'

            acc = xc*yc';

            x(j,fields)=nus1-acc;

            acc2 = xc*ynv';            
            x(j,fields-1)=S(j)-acc2;
        end

        %Last row of X that couldn't be covered in the previous loop
        xc=X(records,:);
        acc2 = xc*ynv';            
        x(records,fields-1)=S(records)-acc2;

        % %boundary condition for the last column of the matrix U
        % for j = 1:records
        % 
        %     xc=X(j,:);
        %     ynv = Y(NV,:);
        % 
        %     %in acc I need to store the cross 
        %     %product between X rows and the last column of Y'
        % 
        %     acc = xc*ynv';            
        %     x(j,fields-1)=S(j)-acc;
        % end
        

        % for i = 1:fields-1
        %     x(1,i)=0;
        % end

        for i = 1:fields
            x(1,i)=0;
            %x(records,i)=x(records-1,i)+(ds/(fields));
        end

        
        %Boundary Conditions end

        sz = size(y);
        %this contains the volatility steps
        records = sz(1);

        %this contains the rank
        fields = sz(2);

        %testing the difference before applying operators
        MU = x*y';
        %%%%%MU = newXBC*newYBC';
        testBeforeOperator = copyOfU-MU;

        ORIGINAL = X*Y';
        testORIGINAL = copyOfU-ORIGINAL;
        test=0;

        %Rewrite "Derivative" to take a matrix instead of vectors
        %The data are expected in columns
        %One field --> One column
        %Rows contain the datapoint for all the fields (Columns)
        fdX = MDerivative(x,ds);     %first derivative of X
        fdY = MDerivative(y,dv);     %first derivative of Y
        sdX = MSecondDerivative(x,ds);  %second derivative of X
        sdY = MSecondDerivative(y,dv);  %second derivative of Y
    
        spVectorSquare = S'.^2;

        %PWProduct needs to take into account that I am
        %now dealing with X and Y as matrixes and not vectors
        s2sdX=PWProduct(spVectorSquare,sdX);
        vY=PWProduct(V,y);
        vsdY=PWProduct(V,sdY);
        vfdY=PWProduct(V,fdY);
        sfdX=PWProduct(S,fdX);
        thetaKappaElement = PWProduct(kappa*(theta-V),fdY);
        
        %firstElement = 0.5*s2sdX*vY';
        X1=0.5*s2sdX;
        Y1=vY;

        %secondElement = 0.5*sigma^2*x*vsdY';
        X2=0.5*sigma^2*x;
        Y2=vsdY;

        % %thirdElement = sigma*sfdX*vfdY';
        X3 = sigma*sfdX; %*0.25; %added 0.25 because excluded from direct formula
        Y3 = vfdY;
        
        %fourthElement = (r-q)*sfdX*Y';
        X4 = (r-q)*sfdX;
        Y4 = y;

        X5 = x;
        Y5 = thetaKappaElement;

        X6 = -r*x;
        Y6 = y;

        dtsrt = sqrt(dt);

        % XC = dtsrt*[X1,X2,X3,X4,X5,X6];
        % YC = dtsrt*[Y1,Y2,Y3,Y4,Y5,Y6];

        XC = dtsrt*[X1,X2,X3,X4,X5,X6];
        YC = dtsrt*[Y1,Y2,Y3,Y4,Y5,Y6];

        newX = [x XC];
        newY = [y YC];

        MU = newX*newY';
        %%%%%MU = newXBC*newYBC';
        testBeforeCompression = U-MU;
        test =1;

        try
            [X,Y]=CompressData(newX,newY,epsilon);
        catch
            a = 1;
        end
        
        %MU = X*Y';
        % test = U-MU;
        %fprintf('Reconstruction Error for the approximated U matrix: %f\n', norm(U-MU)/norm(U));
    end

end

function [X, Y]=CompressData(XC,YC, epsilon)

    [Qx,Rx]=qr(XC,0);
    [Qy,Ry]=qr(YC,0);

    Ri = Rx*Ry';
    
    rank=SVDRank(Ri,epsilon);

    %all the below should be part of SVDRank
    %for learning pupose I keep it in this method
    %once it goes into production, the code will be
    %refactored, consolidated and tested

    [u,s,v]=svd(Ri);

    %the prefix d stands for "dominant"
    du = u(:,1:rank);
    ds = s(1:rank,1:rank);
    %dv = v(1:rank,:);   
    dv = v(:,1:rank);   

    %creation of compressed matrices
    X=Qx*du*sqrt(ds);
    Y=Qy*dv*sqrt(ds);

    %creation of compressed matrices
    %wrong code not working
    %X=Qx*du*ds;
    %Y=Qy*dv';
    
end


%writematrix(matrix_rank, 'G:\My Drive\BO\PhD\Heston\RankAnalysis\matrix_rank.csv')
%writematrix(matrix_rank_home_made, 'G:\My Drive\BO\PhD\Heston\RankAnalysis\matrix_rank_home_made.csv')

%Task: rewrite the code for 1 step only (t=0->t=1)
%In the new code I don't use 2 dimensional U and u: it has to work with
%2 vectors that are 1 dimensional (1 row and 1 column).
%the start is U that is 2 vectors, 1 depending on S and one depending on V.
%Like it was done for the boundaries. The stuff that is currently written
%in the two dimensional case, needs to work on vectors rather than 2
%dimensional matrix. One step implementation will show the linear
%complexity leading to the same result of the quadratic complexity (current
%implementation).


function CP = PWProductWRONG(vector, matrix)
    sz = size(matrix);
    records = sz(1);
    fields = sz(2);

    CP=zeros(records,1);

    for i = 1:records
        acc = 0;
        for j = 1:(fields)
            acc=acc+matrix(i,j);
        end
        CP(i)=acc*vector(i);
    end
end

function CP = PWProduct(vector, matrix)
    sz = size(matrix);
    records = sz(1);
    fields = sz(2);

    CP=zeros(records,fields);

    for i = 1:records
        for j = 1:(fields)
            CP(i,j)=matrix(i,j)*vector(i);
        end
    end
end

function CP = PWProductVector(vector1, vector2)
    sz = size(vector1);
    records = sz(2);

    CP=zeros(records);

    for i = 1:records
        CP(i)=vector2(i)*vector1(i);
    end
end

function CP = SqueezeMatrix(matrix)

    sz = size(matrix);
    records = sz(1);
    fields = sz(2);

    CP=zeros(records,1);

    for i = 1:records
        acc = 0;
        for j = 1:(fields)
            acc=acc+matrix(i,j);
        end
        CP(i)=acc;
    end
end


function D = MDerivative(L,delta)

    sz = size(L);

    records = sz(1);
    fields = sz(2);

    %D=zeros(length(L),1);
    D=zeros(records,fields);

    for i = 2:(records-1)
        for j = 1:(fields)
            D(i,j)=(L(i+1,j)-L(i-1,j))/(delta*2);
        end
    end

    %D(1)=(L(2)-L(1))/delta;
    %D(length(L))=(L(length(L))-L(length(L)-1))/delta;
    
    % for j = 1:(fields)
    %     D(1,j)=(L(2,j)-L(1,j))/delta;
    %     D(records,j)=(L(records,j)-L(records-1,j))/delta;
    % end
end

function D = MSecondDerivative(L,delta)

    sz = size(L);

    records = sz(1);
    fields = sz(2);

    %D=zeros(length(L),1);
    D=zeros(records,fields);

    for i = 2:(records-1)
        for j = 1:(fields)
            D(i,j)=(L(i+1,j)+L(i-1,j)-2*L(i,j))/(delta^2);
        end
    end

    %D(1)=(L(2)-L(1))/delta;
    %D(length(L))=(L(length(L))-L(length(L)-1))/delta;
    
    % for j = 1:(fields)
    %     D(1,j)=(L(2,j)-L(1,j))/delta;
    %     D(records,j)=(L(records,j)-L(records-1,j))/delta;
    % end
end

%derS=(u(s+1,v)-u(s-1,v))/(2*ds);
%derS2=(u(s+1,v)-2*u(s,v)+u(s-1,v))/(ds^2);




function R = SVDRank(U, epsilon)

    %builtin_rank = rank(U,epsilon);

    [uc,sc,vs]=svd(U);

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

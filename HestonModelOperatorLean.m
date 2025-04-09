function [A,B] = HestonModelOperatorLean(A1KX,A1KY,A2KX,A2KY,A3KX,A3KY,A4KX,A4KY,A5KX,A5KY)
    A(:,:,1)=A1KX;
    A(:,:,2)=A2KX;
    A(:,:,3)=A3KX;
    A(:,:,4)=A4KX;
    A(:,:,5)=A5KX;

    B(:,:,1)=A1KY;
    B(:,:,2)=A2KY;
    B(:,:,3)=A3KY;
    B(:,:,4)=A4KY;
    B(:,:,5)=A5KY;

    % A(:,:,1)=sparse(A1KX);
    % A(:,:,2)=sparse(A2KX);
    % A(:,:,3)=sparse(A3KX);
    % A(:,:,4)=sparse(A4KX);
    % A(:,:,5)=sparse(A5KX);
    % 
    % B(:,:,1)=sparse(A1KY);
    % B(:,:,2)=sparse(A2KY);
    % B(:,:,3)=sparse(A3KY);
    % B(:,:,4)=sparse(A4KY);
    % B(:,:,5)=sparse(A5KY);



end

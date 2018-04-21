function [inv] = LUinv(A)
%   LUinv finds the inverse of a square matrix by method of LU
%   decomposition. It depends upon the LUdec() function (see LUdec for
%   more information on LU decomposition.
%
%   Inputs:
%       A           =   square matrix
%
%   Outputs:
%       inv         =   inverse of square matrix A
%
%   Authors:
%       Reed Bell   -   rbel068@aucklanduni.ac.nz

    %Initialisation
    [m,n] = size(A);
    I = eye(n);
    y = zeros(n,m);
    x = y;

    %Call LUdec() to find lower, upper and pivot matrices
    [L, U, P] = LUdec(A);
    
    %Since A*inv(A) = I, the inverse can be determined by solving LUinv(A) = I.
    %Thus, we solve the inverse for A by solving for each column.
    for i = 1:m
        
        %Use forward subsitution to find the values of intermediate vector
        %y where Ly(:,i) = I(:,i) NOTE TO SELF: CURRENT METHOD INVOLVES
        %UNECESSARY CALCULATIONS (MULTIPLICATION OF ENTIRE ROW / COLUMNS).
        for j = 1:m
            y(j,i) = (I(j,i) - (L(j,:)*y(:,i)))/L(j,j);
        end
        
        %Use backward subsitution to find the inverse of A for column i
        %where Ux(:,i) = y(:,i).
        %
        %NOTE TO SELF: INVESTIGATE WHETHER IT IS POSSIBLE TO MERGE THESE
        %LOOPS (CURRENT METHOD RELIES UPON INTERMEDIATE VECTOR BEING
        %CALCULATED FIRST).
        %
        %NOTE TO SELF: IT MIGHT BE CHEAPER TO CALC ALL OF 'y' FIRST THEN
        %CALC ALL OF 'x'.
        for j = m:-1:1
            x(j,i) = (y(j,i) - U(j,:)*x(:,i))/U(j,j);
            
        end
 
    end
    inv = P*x;
end
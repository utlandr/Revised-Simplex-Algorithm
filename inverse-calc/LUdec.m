function [L,U,P] = LUdec(A)
%   LUdec peforms LU decompisition of a square matrix into lower (L) and upper (P)
%   triagular matrices as well as any pivots (P) required such that PA =LU.
%
%   Inputs:
%       A   =   square matrix
%
%   Outputs:
%       L   =   lower triagular matrix 
%       U   =   upper triangular matrix
%       P   =   pivot (permutation) matrix

    %Initialisation
    [m, n] = size(A);
    L = eye(n);
    P = eye(n);
    U = A;

    %Gaussian-Jordan elimination iterations
    for i = 1:m-1
        %Itentify row with the largest magnitude element in ith column
        piv = max(abs(U(i:m,i)));
        
        for k = i:m
            if(abs(U(k,i)) == piv)
                piv_ind = k;
                break;
            end
        end
        %Swap row i with the pivot row in U (in the upper end) 
        U([i,piv_ind], i:m) = U([piv_ind,i], i:m);

        %Same for L (in the lower end)
        L([i,piv_ind], 1:i-1) = L([piv_ind,i], 1:i-1);

        %Same for P (across all rows)
        P([i,piv_ind], :) = P([piv_ind,i], :);

        %Perform elimination
        for j = i+1:m
            L(j,i) = U(j,i)/U(i,i);
            U(j,i:m) = U(j,i:m) - L(j,i)*U(i,i:m);
        end
    end
 
end
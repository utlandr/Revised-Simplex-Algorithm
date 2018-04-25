function [s,minrc] = fullfindEV(n,c,A,varstatus,pi,phase1)
% Returns the index of the entering variable and it's reduced cost,
% or returns 0 if no entering variable exists
% 
%   Input:
%       n = number of variables
%       c = nx1 cost vector
%       A = mxn constraint matrix
%       varstatus = 1xn vector, varstatus(i) = position in basis of variable i, or 0 if variable i is nonbasic
%       pi = mx1 vector of shadow prices
%       phase1 = boolean, phase1 = true if Phase 1, or false otherwise
%
%   Output:
%       s = index of the entering variable
%       minrc = reduced cost of the entering variable
%
%   Author:
%       Reed Bell   -   rbel068@aucklanduni.ac.nz
    
    %Form N from constraint matrix as well as necessary parameters
    nbasicvars = find(varstatus == 0);
    N = A(:,nbasicvars); 
    Nlen = size(N,2);
    
    %Compute vector of reduced costs for non-basic variables. This is
    %looped for proof of understanding. However, the matrix operation below
    %is much more efficient:
    %=======================================================
    %redCost = transpose(c(varstatus==0)) - transpose(pi)*N;
    %=======================================================
    
    %Pre-allocation if in phase1, cN is zero for non-artificial variables
    if phase1
        cN = c(nbasicvars);
        cN(find(nbasicvars) <= n ) = 0;
       
    else
        cN = c(nbasicvars);
        
    end
    
    %Pre-allocation before looping
    m= length(pi);
    piT = transpose(pi);
    redCost = zeros(1,Nlen);
    
    %Iterate through and calculate reduced cost for each variable (and form
    %reduced costs vector
    for i=1:Nlen
        
        %Perform matrix multiplication to get piTN for each column
        piTN = 0;
        for j = 1:m
            piTN = piTN + piT(j)*N(j,i); 
        
        end
        
        %Calculate reduced cost for variable (cN - piT*N)
        redCost(i) = cN(i) - piTN;
        
    end
    
    %Find the index of the variable (x_s) that will reduce the objective function the
    %most (if all are non-negative, then return values [minrc,s] = [0,0]
    [minrc,sPos] = min(redCost);
    s = nbasicvars(sPos);
    
    if minrc >= 0 
        minrc = 0;
        s = 0;
    
    end
    
end
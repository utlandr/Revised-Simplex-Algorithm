function [s, minrc] =   GJfindEV(m,n, c, A, varstatus,pi)
%   Returns the index of the entering variable and it's reduced cost,
%   or returns 0 if no entering variable exists
%
%   Input:
%       m,n         =   number of constraints and variables
%       c           =   nx1 cost vector
%       A           =   mxn constraint matrix
%       varstatus   =   1xn matrix containing position in basis of x_i,
%                       zero if not in the basis
%       pi          =   mx1 dual vector
%
%   Output:
%       s           =   index of the entering variable
%       minrc       =   reduced cost of the entering variable
%
%   Author:
%       Reed Bell   -   rbel068@aucklanduni.ac.nz
    
    %Form N from constraint matrix 
    N = A(:,varstatus == 0); 
    
    %Compute vector of reduced costs for non-basic variables. 
    redCost = transpose(c(varstatus==0)) - transpose(pi)*N;
    
    %Find the index of the variable (x_s) that will reduce the objective function the
    %most. If no variables decrease the objective function, then the current
    %solution is optimal.
    [minrc,s] = min(redCost);
    
    if minrc >= 0 
        minrc = 0;
        s = 0;
    
    end
end
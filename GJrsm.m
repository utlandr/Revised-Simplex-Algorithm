function [result,z,x,pi,basicvars]  =   GJrsm(m,n,c,A,b,basicvars)
%   Solves a linear program using the Revised Simplex Method
%   Assumes standard computational form
%   Starts from the given initial basic feasible solution
%   
%   Inputs:
%       m,n         =   number of constraints and variables
%       c           =   nx1 cost vector
%       A           =   mxn constraint matrix
%       b           =   mx1 rhs vector
%       basicvars   =   1xm vector of indices of basic variables
%       m,n         =   number of constraints and variables
%
%   Outputs:
%       result      =   1 if problem is optimal, -1 if unbounded
%       z           =   objective function value
%       x           =   nx1 solution vector
%       pi          =   mx1 dual vector
%
%   Author:
%       Reed Bell   -   rbel068@aucklanduni.ac.nz


    loopCheck = 1;
    
    %Partition the constraint matrix (A) into basic and non-basic matrices.
    nbasicvars = (1:n);
    nbasicvars(basicvars) = [];
    B = A(:,basicvars);
    
    %Compute the inverse of the Basis matrix (B) using LU decomposition with
    %partial pivoting (see LUinv and its dependency LUdec for more information)
    Binv = LUinv(B);
    
    while loopCheck 
        
        %Compute pi. That is, the vector of duals or shadow prices. 
        pi = transpose(transpose(c(basicvars))*Binv);
        
        %Setup variable status to determine position of each variable in
        %the matrix
        [~,varstatus] = ismember(1:n,basicvars);
        
        %Call findEV() to determine minimum reduced cost
        [s, minrc] = GJfindEV(m,n,c,A,varstatus,pi);
        
        %Condition for optimality
        if s == 0 
            %Notify user and return values
            notify = sprintf("Optimal Found\n");
            loopCheck = 0;
            result = 1;
            x = Binv*b;
            z = transpose(c(basicvars))*x;
            
        else
            BinvAs = Binv*A(:,nbasicvars(s));
            xB = Binv*b;
            [r,minratio] = GJfindLV(m,xB,BinvAs);
            
            if r == 0
                notify = sprintf("Problem is unbounded\n");
                loopCheck = 0;
            
            else
                %Call update() before next iteration to set next iteration
                [varstatus,basicvars, cB, Binv, xB]  =  GJupdate(m, c, s, r, BinvAs, varstatus, basicvars, c(basicvars), Binv, xB);             
            
            end
        end
    end
    fprintf(notify) 
end                
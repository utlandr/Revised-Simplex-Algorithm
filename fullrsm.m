function [result,z,x,pi]  =   fullrsm(m,n,c,A,b)
%   Solves a linear program using the Revised Simplex Method
%   Assumes standard computational form
%   Starts from the given initial basic feasible solution
%   
%   Inputs:
%       m,n         =   number of constraints and variables
%       c           =   nx1 cost vector
%       A           =   mxn constraint matrix
%       b           =   mx1 rhs vector
%       m,n         =   number of constraints and variables
%
%   Outputs:
%       result      =   1 if problem is optimal, -1 if unbounded
%       z           =   objective function value
%       x           =   nx1 solution vector
%
%   Author:
%       Reed Bell   -   rbel068@aucklanduni.ac.nz

    
    %   INITIALISATION STAGE
    %   =====================
    
    loopCheck = 1;
    
    %Partition the constraint matrix (A) into basic matrix.
    nbasicvars = (1:n);
    nbasicvars(basicvars) = [];
    B = A(:,basicvars);
    
    %Setup variable status to determine position of each variable in
    %the matrix
    [~,varstatus] = ismember(1:n,basicvars);
    
    %Compute the inverse of the Basis matrix (B) using LU decomposition with
    %partial pivoting (see LUinv and its dependency LUdec for more information)
    Binv = LUinv(B);
    
    %Compute basic variable values
    xB = Binv*b;
    
    %Get basic coefficients subset
    cBT = transpose(c(basicvars));
    
    
    %   RSM ITERATION
    %   =============
    
    while loopCheck 
        
        %Compute pi. That is, the vector of duals or shadow prices. 
        pi = transpose(cBT*Binv);
               
        %Call GJfindEV() to determine minimum reduced cost
        [s, ~] = GJfindEV(m,n,c,A,varstatus,pi);
        
        %Condition for optimality
        if s == 0 
            %Notify user and return values
            notify = sprintf("Optimal Found\n");
            loopCheck = 0;
            result = 1;
            x = Binv*b;
            z = cBT*x;
            
        else
            BinvAs = Binv*A(:,nbasicvars(s));
            
            %call GJfindLV() to determine leaving variable
            [r,~] = GJfindLV(m,xB,BinvAs);
            
            if r == 0
                notify = sprintf("Problem is unbounded\n");
                loopCheck = 0;
            
            else
                %Call GJupdate() before next iteration to set next iteration
                [varstatus,basicvars, cB, Binv, xB]  =  GJupdate(m, c, s, r, BinvAs, varstatus, basicvars, transpose(cBT), Binv, xB);             
            
            end
        end
    end
    fprintf(notify) 
end                
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
    
    %With no initial basis provided, bootstrapping is first required to
    %create a basis before solving the LP.
    phase1 = true;
    
    %Identify basic and non-basic variables. Because this is
    %initialisation, all basic variables will be artificial (implicity
    %implied by 1:n+m rather than 1:n. Basic and non-basc constraint 
    %matrices can be setup(as they are grouped together).
    basicvars = n+1:n+m;
    nbasicvars = 1:n;
    B = eye(m);
    N = A; 
    A = [N,B];
    
    %Setup variable status to determine position of each variable in
    %the matrix
    [~,varstatus] = ismember(1:n+m,basicvars);
    
    %Compute the inverse of the Basis matrix (B). The inverse is the 
    %identity matrix (all variables are artificial thus B = I) 
    Binv = B;
    
    %Compute basic variable values
    xB = Binv*b;
    
    %Get basic coefficients subset
    c = [c;ones(m,1)];
    cBT = c;
    cBT = transpose(cBT(basicvars)); 
    
    
    %   RSM ITERATION
    %   =============
    
    %Begin RSM iterations to minimise objective function. If phase1 is true
    %then we begin phase 1 to create a viable basis.
    while loopCheck 
        
        %Compute pi. That is, the vector of duals or shadow prices. 
        pi = transpose(cBT*Binv);
               
        %Call GJfindEV() to determine minimum reduced cost
        [s, ~] = fullfindEV(n,c,A,varstatus,pi,phase1);
        
        %Condition for optimality
        if s == 0 
            %return values
            %notify = sprintf("Optimal Found\n");
            loopCheck = 0;
            result = 1;
            x = Binv*b;
            z = cBT*x;
            
        else
            BinvAs = Binv*A(:,nbasicvars(s));
            
            %call GJfindLV() to determine leaving variable
            [r,~] = GJfindLV(m,xB,BinvAs, phase1, basicvars);
            
            if r == 0
                %notify = sprintf("Problem is unbounded\n");
                loopCheck = 0;
            
            else
                %Call GJupdate() before next iteration to set next iteration
                [varstatus,basicvars, cB, Binv, xB]  =  GJupdate(m, c, s, r, BinvAs, varstatus, basicvars, transpose(cBT), Binv, xB);             
            
            end
        end
    end
    fprintf(notify) 
end                
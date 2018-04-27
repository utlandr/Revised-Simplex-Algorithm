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
%       x           =   nx1 solution vector arranged such that x(i) is x_i
%       pi          =   shadow costs of the basic variables. 
%
%   Author:
%       Reed Bell   -   rbel068@aucklanduni.ac.nz


    %   INITIALISATION STAGE
    %   =====================
    
    %looping condition (
    result = 2;
    
    %With no initial basis provided, bootstrapping is first required to
    %create a basis before solving the LP.
    phase1 = true;
    
    %Identify basic and non-basic variables. Because this is
    %initialisation, all basic variables will be artificial (implicity
    %implied by 1:n+m rather than 1:n. Basic and non-basc constraint 
    %matrices can be setup(as they are grouped together).
    basicvars = n+1:n+m; B = eye(m); N = A; A = [N,B];
    
    %Setup variable status to determine position of each variable in
    %the matrix
    [~,varstatus] = ismember(1:n+m,basicvars);
    
    %Compute the inverse of the Basis matrix (B). The inverse is the 
    %identity matrix (all variables are artificial thus B = I) 
    Binv = B;
    
    %Compute basic variable values
    xB = Binv*b;
    
    %Reform c for phase 1 and create 
    cB = ones(m,1);
    c = [c;cB];
    cBT = transpose(cB); 
    %searchC = find( varstatus ~= 0);
    %cBT(varstatus(searchC( searchC<= n))) = 0;    
    
    %   RSM ITERATION
    %   =============
    
    %Begin RSM iterations to minimise objective function. If phase1 is true
    %then we begin phase 1 to create a viable basis.
    while result == 2 
        
        %Compute pi. That is, the vector of duals or shadow prices. 
        pi = transpose(cBT*Binv);
        
        %Call fullfindEV() to determine minimum reduced cost
        [s, ~] = fullfindEV(n,c,A,varstatus,pi,phase1);
        
        %Condition for optimality (empty set when A matrix is square)
        if isequal(s, 0) || isempty(s)
            %return values
            if ~phase1
                result = 1;
            end
            
            %Calculate solution and objective values (whilst also
            %rearranging x such that x(i) is the value for x_i
            x = zeros(n,1);
            x(varstatus ~= 0) = xB(varstatus(varstatus ~= 0));
            z = cBT*xB;    
            
            
            if phase1
                %Check for positive objective (infeasible)
                if z > 0
                    result = 0;
                    z = NaN;
                    x = [];
                    pi = [];
                    
            break;
            
                else
                    %End of phase 1
                    phase1 = false;
                    
                    %Remove all artificial variables not in the basis
                    nsearch = find(varstatus ==0);
                    varstatus(nsearch(nsearch > n)) = [];
                
                end
                
            end
            
        else
            %Calculate BinvA for the entering variable
            BinvAs = Binv*A(:,s);
            
            %Call fullfindLV() to determine leaving variable
            [r,~] = fullfindLV(n,xB,BinvAs, phase1, basicvars);
            
            %Check if fullfindLV discovers unboundedness 
            if r == 0
                result = -1;
                z = NaN;
                pi = [];
                x = [];
                break;
            else
                %Call fullupdate() before next iteration to set next iteration
                [varstatus,basicvars, cB, Binv, xB]  =  fullupdate(m, c, s, r, BinvAs, phase1, varstatus, basicvars, transpose(cBT), Binv, xB);
                
                %End phase 1 if there are no artificial variables in the
                %basis (or are degenerate)
                if phase1 && (transpose(cB)*xB == 0)
                   phase1 = false;
                   cB = zeros(m,1); cB(varstatus(varstatus>0)) = c(varstatus>0);
                end
                
                cBT = transpose(cB);    
            end
            
        end
        
    end
     
end                
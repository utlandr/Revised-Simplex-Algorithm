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
    N = A(:,nbasicvars);
    
    %Compute the inverse of the Basis matrix (B) using LU decomposition with
    %partial pivoting (see LUinv and its dependency LUdec for more information)
    Binv = LUinv(B);
    
    while loopCheck 
        %Initialisation for ratio test
        ratioTest = zeros(2,m);

        %Compute pi. That is, the vector of duals or shadow prices. 
        tranPi = transpose(c(basicvars))*Binv;

        %Compute vector of reduced costs for non-basic variables. 
        redCost = transpose(c(nbasicvars)) - tranPi*N;

        %Find the index of the variable (x_s) that will reduce the objective function the
        %most. If no variables decrease the objective function, then the current
        %solution is optimal.
        [val,s] = min(redCost);

        if val >= 0 
            %Notify user and return values
            notify = sprintf("Optimal Found\n");
            loopCheck = 0;
            
            result = 1;
            x = Binv*b;
            z = transpose(c(basicvars))*x;
            pi = transpose(tranPi);
        else

            %Apply minimum ratio to determine the leaving variable index (r).
            ratioTest(1,:) = transpose(Binv*b);
            ratioTest(2,:) = Binv*A(:,nbasicvars(s));
            oldVars = transpose(ratioTest);

            %If all values for entering variable are <=0, the problem is
            %unbounded.
            if ratioTest(2,:) <= 0
                notify = sprintf("Problem is unbounded\n");
                loopCheck = 0;
            else

            %Remove values for the entering variable that are not positive
            ratioTest(:,ratioTest(2,:)<=0) = NaN;
            [~,r]   = min(ratioTest(1,:)./ratioTest(2,:));


            %Perform Gaussian-Jordan Pivoting replace variables and find new Binv
            augMatrixOld = [oldVars(:,1), Binv, oldVars(:,2)];

            %Normalise row reprsenting leaving variable
            normrMatrix  = augMatrixOld;
            normrMatrix(r,:)  = normrMatrix(r,:)/normrMatrix(r,end);

            elimRows = 1:m;
            elimRows(elimRows==r) = [];
            for i = elimRows

                normrMatrix(i,:) = normrMatrix(i,:) - (normrMatrix(i,end)/normrMatrix(r,end))*(normrMatrix(r,:));

            end

            %Update vector tracking basic and non basic variables (in order)
            tmp = basicvars;
            basicvars(s) = nbasicvars(r);
            nbasicvars(r) = tmp(s);

            %Update B and N 
            B = A(:,basicvars);
            N = A(:,nbasicvars);
            
            Binv = normrMatrix(:,2:end-1);

            end
        end
    end
    fprintf(notify)
end



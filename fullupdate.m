function [varstatus,basicvars,cB,Binv,xB] = fullupdate(m,c,s,r,BinvAs,phase1,varstatus,basicvars,cB,Binv,xB)
% Updates the basis representation.
%
%   Input:
%       m = number of constraints
%       c = nx1 cost vector
%       s = index of entering variable
%       r = position in the basis of the leaving variable
%       BinvAs = mx1 Binv*As vector
%       phase1 = boolean, phase1 = true if Phase 1, or false otherwise
%       varstatus = 1xn vector, varstatus(i) = position in basis of variable i,
%                   or 0 if variable i is nonbasic
%       basicvars = 1xm vector of indices of basic variables
%       cB = mx1 basic cost vector
%       Binv = mxm basis inverse matrix
%       xB = mx1 basic variable vector
%
%   Output:
%       varstatus = 1xn updated varstatus vector
%       basicvars = 1xm updated basicvars vector
%       cB = mx1 updated basic cost vector
%       Binv = mxm updated basis inverse matrix
%       xB = mx1 updated basic variable vector%
%
%   Author:
%       Reed Bell   -   rbel068@aucklanduni.ac.nz

    %Setup augmented matrix
    nbasicvars = find(varstatus==0);
    augMatrix = [xB, Binv, BinvAs];
    augMatrix(r,:)  = augMatrix(r,:)/augMatrix(r,end);
    
    %Identify rows to pivot against
    elimRows = 1:m;
    elimRows(elimRows==r) = [];
    
    %Perform Gaussian Jordan pivoting to obtain new inverse basis and x solution. 
    for i = elimRows
        augMatrix(i,:) = augMatrix(i,:) - (augMatrix(i,end)/augMatrix(r,end))*(augMatrix(r,:));

    end

    %Update vector tracking basic and non basic variables (in order).
    if basicvars(r) > (length(c) - m)
        artVar = basicvars(r);
        [basicvars(r),nbasicvars(nbasicvars == s)] = deal(s,basicvars(r));
        [varstatus(artVar),varstatus(s)] = deal(0,r) ;
        
        
    else
        [varstatus(basicvars(r)),varstatus(s)] = deal(0,r);
        [basicvars(r),nbasicvars(nbasicvars == s)] = deal(s,basicvars(r));
    
    end
    
    %Check for completion of phase 1 (no artificial variables in
    %the basis).
    n = length(varstatus) - m;
    if phase1 && isempty(basicvars(basicvars > n))
        %Removal of artificial variables that are non basic
        varstatus(nbasicvars(nbasicvars > n)) = [];
        
    end

    %Update function return values (phase 1 non-artificial variables
    %set to 0).
    cB(r) = c(s)*~phase1;
    Binv = augMatrix(:,2:end-1);
    xB = augMatrix(:,1);
    
end
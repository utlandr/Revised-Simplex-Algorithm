function [varstatus,basicvars,cB,Binv,xB]  =  GJupdate(m,c,s,r,BinvAs,varstatus,basicvars,cB,Binv,xB)
%   Updates the basis representation.
%
%   Input:
%       m,n         =   number of constraints and variables
%       c           =   nx1 cost vector
%       A           =   mxn constraint matrix
%       b           =   mx1 rhs matrix
%       s           =   index of entering variable
%       r           =   position in the basis of the leaving variable
%       basicvars   =   1xm vector of indices of basic variables
%       cB          =   mx1 basic cost vector
%       B           =   mxm basis matrix
%       xB          =   mx1 basic variable vector
%
%   Output:
%       basicvars   =   1xm updated basicvars vector
%       cB          =   mx1 updated basic cost vector
%       B           =   mxm updated basis matrix
%       xB          =   mx1 updated basic variable vector
%
%   Author:
%       Reed Bell   -   rbel068@aucklanduni.ac.nz

    %Perform Gaussian-Jordan Pivoting replace variables and find new Binv
    nbasicvars = find(varstatus==0);
    augMatrix = [xB, Binv, BinvAs];

    %Normalise row representing leaving variable
    augMatrix(r,:)  = augMatrix(r,:)/augMatrix(r,end);

    elimRows = 1:m;
    elimRows(elimRows==r) = [];
    for i = elimRows
        augMatrix(i,:) = augMatrix(i,:) - (augMatrix(i,end)/augMatrix(r,end))*(augMatrix(r,:));

    end

    %Update vector tracking basic and non basic variables (in order)
    basicvars(s) = nbasicvars(r);
    [~,varstatus] = ismember(1:size(c),basicvars);
    
    %Update function return values
    cB = c(basicvars);
    Binv = augMatrix(:,2:end-1);
    xB = augMatrix(:,1);
    
end
function [r, minratio]  =  fullfindLV(m, xB, BinvAs, phase1, basicvars)
%   Returns the position in the basis of the leaving variable,
%   or returns 0 if no leaving variable exists
% 
%   Input:
%       m           =   number of constraints and variables
%       xB          =   mx1 basic variable vector
%       BinvAs      =   mx1 vector of Binv*As
%       phase1      =   boolean, phase1 = true if Phase 1, or false otherwise
%       basicvars   =   1xm vector of indices of basic variables
%
%   Output:
%       r           =   position in the basis of the leaving variable
%       minratio    =   minimum ratio from ratio test
%
%   Author:
%       Reed Bell   -   rbel068@aaucklanduni.ac.nz
    
    %Setup ratio test matrix
    ratioTest = [xB,BinvAs];

    %If all values for entering variable are <=0, the problem is unbounded
    if ratioTest(:,2) <= 0
        r = 0;
        minratio = 0; 
    
    else
        %Remove values for the entering variable that are not positive
        ratioTest(ratioTest(:,2)<=0,:) = NaN;
        [minratio,r]   = min(ratioTest(:,1)./ratioTest(:,2));
    
    end
    
end
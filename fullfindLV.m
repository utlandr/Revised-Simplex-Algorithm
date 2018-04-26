function [r, minratio]  =  fullfindLV(n, xB, BinvAs, phase1, basicvars)
%   Returns the position in the basis of the leaving variable,
%   or returns 0 if no leaving variable exists
% 
%   Input:
%       n           =   number of variables
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
%       Reed Bell   -   rbel068@aucklanduni.ac.nz
    
    %If all values for entering variable are <=0, the problem is unbounded
    if BinvAs <= 0
        r = 0;
        minratio = 0; 
    
    else
        %Setup checks for non-zero artificial variables in the basis
        artCheck = (find(basicvars > n));
        artCheck = find(xB(artCheck) ~= 0);
        
        %Take the first artificial variable when in phase 2
        if ~phase1 && ~isempty(artCheck)
            r = artCheck(0);
            minratio = 0;
    
        else 
            %Setup ratio test matrix
            ratioTest = [xB,BinvAs];
    
            %Remove values for the entering variable that are not positive
            ratioTest(ratioTest(:,2)<=0,:) = NaN;
            [minratio,r]   = min(ratioTest(:,1)./ratioTest(:,2));
      
        end
        
    end
    
end
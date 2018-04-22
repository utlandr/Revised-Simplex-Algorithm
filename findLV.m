function [r, minratio]  =  findLV(m, xB, BinvAs)
%   Returns the position in the basis of the leaving variable,
%   or returns 0 if no leaving variable exists
% 
%   Input:
%       m           =   number of constraints and variables
%       xB          =   mx1 basic variable vector
%       BinvAs      =   mx1 vector of Binv*As
%
%   Output:
%       r           =   position in the basis of the leaving variable
%       minratio    =   minimum ratio from ratio test
%
%   Author:
%       Reed Bell   -   rbel068@aaucklanduni.ac.nz
    
    %Initialisation for ratio test
    ratioTest = zeros(2,m);

    %Apply minimum ratio to determine the leaving variable index (r).
    ratioTest(1,:) = transpose(Binv*b);
    ratioTest(2,:) = transpose(BinAs);%Binv*A(:,nbasicvars(s));
    oldVars = transpose(ratioTest);

    %If all values for entering variable are <=0, the problem is
    if ratioTest(2,:) <= 0
        r = 0;
        minratio = 0; 
    
    else
        %Remove values for the entering variable that are not positive
        ratioTest(:,ratioTest(2,:)<=0) = NaN;
        [minratio,r]   = min(ratioTest(1,:)./ratioTest(2,:));
    
    end
    
end
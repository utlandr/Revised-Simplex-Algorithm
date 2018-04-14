function [r, minration]  =  findLV(m, xB, BinvAs)
%   Returns the position in the basis of the leaving variable,
%   or returns 0 if no leaving variable exists
% 
%   Input:
%       m           =   umber of constraints and variables
%       xB          =   mx1 basic variable vector
%       BinvAs      =   mx1 vector of Binv*As
%
%   Output:
%       r           =   position in the basis of the leaving variable
%       minratio    =   minimum ratio from ratio test
%
%   Author:
%       Reed Bell   -   rbel068@aaucklanduni.ac.nz
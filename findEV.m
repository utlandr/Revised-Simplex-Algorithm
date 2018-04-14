function [s, minrc] =   findEV(n, c, A, pi)
%   Returns the index of the entering variable and it's reduced cost,
%   or returns 0 if no entering variable exists
%
%   Input:
%       m,n         =   number of constraints and variables
%       c           =   nx1 cost vector
%       A           =   mxn constraint matrix
%       pi          =   mx1 dual vector
%
%   Output:
%       s           =   index of the entering variable
%       minrc       =   reduced cost of the entering variable
%
%   Author:
%       Reed Bell   -   rbel068@aucklanduni.ac.nz
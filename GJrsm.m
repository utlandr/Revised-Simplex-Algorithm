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
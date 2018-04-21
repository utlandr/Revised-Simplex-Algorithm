function    []  =   TestWrap(mat)
%   TestWrap.m is a wrapper script that tests various stages of the Revised
%   Simplex algorithm implementation using preworked problems. It parses a
%   .mat file containing the necessary data for inputting into GJrsm().
%
%   Inputs:
%       mat         =       .mat files containing necessay data for test
%                           conditions.
%
%   Outputs:
%       
%   Author:
%       Reed Bell   -       rbel068@aucklanduni.ac.nz

%Load .mat file containing test cases
load(mat)

%Perform the Revised Simplex Method

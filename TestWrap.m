function    []  =   TestWrap(testDir,rsmDir,testTypes)
%   TestWrap.m is a wrapper script that tests various stages of the Revised
%   Simplex algorithm implementation using preworked problems. Assumes
%   Assignment1 naming scheme for 1c
%
%   Inputs:
%       testDir         =       Directory string (can be relative or absolute)
%                               containing test case (.m) files. Default is
%                               './test-cases/'
%
%       rsmDir          =       Directory string to RSM algorithm wrapper
%                               script 'GJrsm or fullrsm'. Default is CWD.
%
%       testTypes       =       vector of strings specifying specific tests
%                               or test files to be used. 
%
%   Outputs:
%       
%   Author:
%       Reed Bell   -       rbel068@aucklanduni.ac.nz

%Default input values would be nice, but can't be stuffed atm to user
%matlabs' parser (what appears to be much uglier than python's argparser).
%So, for now, we deal with putting the right value in

%Set working directory to where the RSM implementation is located
cd(rsmDir)

%Grab required .mat files

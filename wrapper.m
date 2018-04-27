%Wrapper script that iterates through folder containing test files and
%performs the RSM algorithm on it.
clear
testFiles = dir('**/test-cases/test*');
iter = struct2cell(testFiles);
iter = iter([1,2,],:);

%Iterate through test cases 
for iter = iter
    
    %Clear space of function input variables
    clearvars -except testFiles iter
    load(char(strcat(iter(2),'/',iter(1))))
    disp(iter(1))
    
    %Call function (trying)
    try
        [result,z,x,pi]  =   fullrsm(m,n,c,A,b)
        
    catch ME
        disp(ME.identifier)
        
    end
    
   
    
end

clear
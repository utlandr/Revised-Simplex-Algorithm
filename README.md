# Matlab Simplex Algorithm

This is an implementation of the Simplex Algorithm in Matlab as part of my studies in Linear Programming. This implementation is written in Matlab and is currently capable of solving tested LP problems whilst using bootsterapping methods to create a basis. Code for LU decomposition and matrix inverse calculations (now redundant) is also included.

###### Disclaimer
This is an educational implementation of the RSM algorithm. Certain edge cases may not have been accounted for. 


## How does it work?
A linear program can described by a 'basis' subset of the full problem that has the same number of elements as there are constraints. These solutions are called 'basic feasible solutions' but are not always optimal.

The Simplex methods iteratively tests different combinations of basic feasible solutions until it discovers conditions only found when optimal. Rather than a brute force method testing every single combination, the simplex method can, in each iteration choose combinations that only improve the solution, or detect conditions that may render the problem unbounded or infeasible.

### Bootstrapping
This implementation of the RSM algorithm includes a bootstrapping phase. The RSM algorithm technically works with the inverse of the basis matrix. One way to find an inverse basis is to 'randomly' select variables that will make the basic feasible solution and then solve the inverse (using LU decomposition or Gaussian Elimination). However, for large matrices, this is computationally expensive and impractical. 

The alternative method is to create artificial variables and use RSM to iteratively drive out each artificial variable (which are replaced with variables in the actual linear program). This creates a basis matrix to then solve for the linear program of interest. 

## Notes

I have also included scripts for finding the inverse of a matrix using LU decomposition. This was before I implemented the bootstrapping method.

There are much better algorithms and solvers available for practical applications. AMPL is capable of using more modern solvers such as CPLEX and GUROBI to solve more complex (and integer) problems. Matlab has its own linear program solver `linprog` that uses the Dual Simplex Method (DSM). 


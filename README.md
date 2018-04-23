# Matlab Simplex Algorithm

This is an implementation of the Simplex Algorithm in Matlab as part of my course studies. This implementation is written in Matlab and is currently capable of solving linear programs with provided basis matrix. The next step is to apply bootstrapping methods to create our own basis matrix.

## How does it work?
A linear program comes in the form:

$$$$

It can be shown that for any linear program in this format, a solution for x can be described by a basis subset that has the same number of elements as there are constraints. These solutions are called 'basic feasible solutions' but are not always optimal.

The Simplex methods iteratively tests different combinations of basic feasible solutions until it discovies conditions only found when optimal. Rather than a brute force method testing every single combination, the simplex method can, in each iteration choose combinations that only improve the solution. 

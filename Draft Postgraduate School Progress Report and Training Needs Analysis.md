Draft Summary in point form of research completed at 1 January 2015
Fixed effects zero-inflated Poisson model
Mixed effects model random intercepts with Inverse Gamma prior
What I did for the next six months was Inverse Wishart
Random intercept
Random slopes
Splines
Four different algorithms that solve the problem
Extensive investigation of performance of various parameterisations of the problem
Solving triangular systems versus inverting entire matrix
Reference BLAS versus Matrix library/sparse matrices versus Eigen/C++. Sparse matrices don't 
help.
We tried many things which didn't work or didn't end up helping:
- custom solvers
- sparse matrices
Fairly deep investigation of R performance optimisation
Parameterisation of Cholesky matrix which ensures postive definiteness
How to speed up Stan sampler so that performance is acceptable - it turns out that of
the equivalent parameterisations, using Cholesky factors is best

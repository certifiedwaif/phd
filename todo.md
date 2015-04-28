To Do list
==========
27/04/2015
----------
- Install R packages
To install limma:
source("http://bioconductor.org/biocLite.R")
biocLite("limma")

To install rstan:
source('http://mc-stan.org/rstan/install.R', echo = TRUE, max.deparse.length = 2000)
install_rstan()

To install the other packages from CRAN:
install.packages(c("Rcpp", "RcppEigen", "numDeriv", "digest"))

23/04/2015
----------
- OH&S - chair
- Boxes for moving
- SU Maker Club - fill out form
- Teach tutorial
	Done
- Check data generation process
	This looks fine. We're sorting the x values appropriately.
- Check MCMC
- Check details of VB
- Meet with John to discuss the problem
	See the meeting notes for 2015 04 23 for more details, but the problem was basically that
	the ZOSull function returns a matrix using the bspline function which includes the basis
	terms [1 x x^2 x^3 (x-K)_+^3 ...]. The first two components [1 x] are exactly the contents
	of $\mX$, so when we form $\mC$, the matrix is not of full rank and so the cross-product
	of mC with itself is not invertible. This of course makes model fitting impossible. It's
	surprising that anything was working at all!
- Clean up mail spool
	Done
- Enter rolls
	Done
- Claim pay
	Done
- Pay phone bill
	Done
	
22/04/2015
----------
- Re-ran MCMC with 20k samples. MCMC and VB still get different \vmu, which can't
	be right.
- Teach tutorial
	Done
- Learned more matrix calculus
	I didn't understand why the adjugate of a matrix behaved as it did. So I proved that
	\mA adj(\mA) = det(A) \mI

21/04/2015
----------
- Figure out why optimiser goes outside of the parameter space.
	Because you didn't tell it not to! Set a lower constraint of -15 for
	everything.
- The VB algorithm stops pretty early. Check that the lower bound is
	really right.
- Check matrix derivative in GVA.

20/04/2015
----------
- Type up notes taken in the meeting this morning with John and Sarah.
	Done, but they need to be checked.
- Continue trying to fit splines.
	The MCMC took forty eight hours, but it completed successfully. Some problems in
	zero_inflated_model.R with not dealing with spline_dim.
	Optimiser is failing with non-finite value.
	The test case I've chosen is pretty terrible, so I should use something simpler.

16/4/2015
---------
- Teach classes, enter rolls
	Done
- The Stan sampler is currently rejecting every proposed sample
	The solution was to have two seperate Stan source files, one for random intercepts and slopes
	and another for splines. There's currently no requirement to have both together, although
	there's no reason why we couldn't do this. It's just harder.
- Read the chapter on non-linear optimisation, particularly as it relates to
	L-BFGS
- Make appointment to see better chair etc.
	I went looking for the relevant room in the IT building, but couldn't find it. I'll go again
	next week.
- Write sigma2 accuracy code
- Read Challis
- Literature review
- Put code into a package. Hadley Wickham shows you how:
	http://r-pkgs.had.co.nz/r.html

15/4/2015
---------
- Continue with splines. Implement something today. A lot of progress was
	made on this
- Read the chapter on non-linear optimisation, particularly as it relates to
	L-BFGS
- Make appointment to see better chair etc.
- Write sigma2 accuracy code
- Read Challis
- Literature review

14/4/2015
---------
- Write up meeting notes from yesterday
- Splines
	Read the relevant sections of Numerical Recipes, my numerical analysis book and 
	semiparametric regression. I have a good understanding now. This ties in well with
	the Computational Projects in Applied Mathematics course that I'm sitting in on.
- Read more Numerical Recipes - more of the linear algebra chapter, and the
	chapter on non-linear optimisation, particularly as it relates to L-BFGS
- Prepare for tutorials
- Implement caching in gaussian.R
- Make appointment to see better chair etc.
- Write sigma2 accuracy code

13/4/2015
---------
- It turns out that Reference BLAS already checks for zeros in its' linear algebra routines.
	So there's no need to optimise this case.
- Present what you've learned to the rest of the group.

9/4/2015
--------
- At the end of yesterday, sparse matrices seemed to finally be helping. Check that this gain is 
	real, not wishful thinking.

	It turns out to be faster than the old algorithm, but only by 10% or so. If you don't use sparse
	matrices at all, the speed of the algorithm immediately doubles. So for small problems, that
	seems to be the way to go. I'd like to try both approaches on a larger problem and see what happens then.

- Continue modifying the rest of the GVA2 code to use the same strategy, and see if the 
	improvements continue.

	I'm afraid that they didn't. In fact, things got worse. This seems to be largely because of the
	S4 machinery you need to go through in order to get to the sparse matrix routines in the Matrix 
	library. I've tried various ways of mitigating this, and nothing's worked thus far.

- Try removing sparse matrices and see what happens then.
	The speed of the code doubles over the old GVA algorithm that we're competing against.

- Keep reading Numerical Recipes to improve understanding. Read numerical linear algebra and 
	optimisation books if there's enough time.


8/4/2015
--------
- Sparse matrices seem to make everything slower. Try to find out why, or if there's a way of 
	speeding them up.
	Ideas to try:
	- different types of problem e.g. different dimension of random effects
	- different ways of constructing the sparse matrices i.e. take advantage of the block structure

7/4/2015
--------
- Literature review - twenty papers on zero-inflated mixed models
- Continue with sparse matrices and solvers in GVA mLambda inverse parameterisation
- Go through gaussian.R eliminating unused parameters, especially mLambda
- Test that derivatives are correct, and that the code really works
- Check speed of new version using a profiler. Perhaps give lineprof a try.
- RLint
- Splines
- Application
- Find a new place to live
- Implement survival analysis model
- Prove fact about p(b_i)
- Mail Ben and SU Maker people about Keiran's idea re: maker activities
- GVA is ~3.2 seconds, whereas GVA2 is currently ~3.7 seconds. Find out where the GVA2 algorithm
  is spending its time.

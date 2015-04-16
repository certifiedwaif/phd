To Do list
==========

14/4/2015
---------
- Write up meeting notes from yesterday
- Splines
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

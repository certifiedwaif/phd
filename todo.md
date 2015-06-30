To Do list
==========

30/06/2015
----------
- Found major bug:
    mSigma.u.inv <- kronecker(diag(1, m - 1), solve(mPsi)/(v - d - 1))
  should have been
    mSigma.u.inv <- kronecker(diag(1, m - 1), solve(mPsi/(v - d - 1)))
  That is, I was inverting and then dividing, rather than dividing and then inverting.

  This had the effect of raising the accuracy of _all_ the intercepts above 80%!

29/06/2015
----------
- Progress update
- Met with John quickly
  John commented that variance component accuracy is low
- Increase degrees of freedom for the Inverse Wishart from 2 to 4. v > p + 1 for the expectation
  to be valid.

28/06/2015
----------
- Install lintr with:
  library(devtools)
  install_github('jimhester/lintr') 
  I haven't got it working yet.
- Bug fixes:
  * 76 GVA2 - Looks like early convergence
    But why? L-BFGS-B can't seem to solve the problem.

25/06/2015
----------
- Using the median accuracy results to find bugs/weird corner cases
  In a couple of cases, optimiser was stopping too early, giving poor results.
  In one case I've found, the optimiser is just going off in a completely
  weird direction and never recovering until finally a non-finite value is return from the
  fn.
  0.5 zeroes is too hard, but we can solve 0.4 zeroes aka rho = 0.6.
- Testing is taking a long time. We should create a set of test files, store them on verona
  or a machine with large disks and then test against those repeatedly. It's currently taking
  ~40 minutes to get feedback, which is stupid.
  Done, and it helped a lot.
- Now I have a lot

23/06/2015
----------
- Work through first four weeks of solutions for Dr Di Warren
  Done.
- Run more median accuracy tasks on Maths department servers
  I'm only running 10,000 MCMC samples for each comparison, which is not very accurate,
  and it's taking ~1/2 a day or more to produce the data.
- Help Jean with parallel R question
  Done. mclapply
- Annotate accuracy plots
  Code written. Testing now.
- accuracy.R should be command-line driven as well?
  Done
- Check when maximum iterations are reached
- Update John
- Fill in and hand in tutoring preferences form
  Holidays in October
- Claim pay for tutoring
- Invoice for MeanPath
- Add legend to accuracy plots
- qqplots?
- Questions: Why does the optimiser go crazy sometimes
  How to check if iterations are exhausted - return the number of iterations used from
  zipvb

22/06/2015
----------
- Finish presentation of tools
- Melanoma meeting
  Jean suggested that I add qqplots, present numerical accuracies
- Give presentation
  That was nerve-wracking
- Meet with John
  Like me, he's nervous about the apparent contradiction between the good accuracy
  plots
- Meeting with John
  Intercept 2
  Slope 1
  rho 0.5
  Variance components .5^2
  Increasing variance components e.g. variance components 3.0 would ruin everything
  Example where accuracy low (i.e. reproducible example, find out why)
  Give data to John for median accuracy
  Things to check:
  * Accuracy poor for all or some
  * GVA NR, is max. iteration limit reached
  * Label plots with accuracy, multiple simulations
  * Same data being supplied to MCMC and VB?
  * Terminating early?
- Found problems:
  * I was initialising vmu with lm or glm in the random slopes and spline cases,
    but not the random intercepts case. This initial value has a large bearing on how
    well the optimisation goes from that point. In median_accuracy, I was using the random
    intercepts case.
  * Data for the random intercepts case was being generated uniformly in the interval
    [-1, 1]. This gives a very narrow range of data, leading to a problem which is very
    hard to solve given how many things we're trying to estimate. For random slopes and splines, was generated from N(0, 1). We now use this option for the random
    intercepts case as well.
  * Low effective sampling of MCMC for some parameters
- Implement changing parameters idea? i.e. stop if lower bound stops increasing and
  parameters stop changing?

21/06/2015
----------
- Write presentation of tools

18/06/2015
----------
- Get accuracy of variance components.
  This turned out to be harder than I would have expected. dgamma() sometimes doesn't
  integrate to 1, and density functions produced by density() sometimes don't either!
- Try to answer more questions from a MATH1015 student, whose exam is tomorrow.

17/06/2015
----------
- Helped a student out with a confidence interval question, because I'm too nice to watch
  them flounder.

16/06/2015
----------
- Rewrote boxplot code for median_accuracy
- Running median accuracy tasks for gva, laplace and gva_nr.
- The accuracy for the slope co-efficient is good, but for all other parameters it's not
  especially good. Perhaps I need to run the MCMC for longer than 1e4.
- coverage_percentage is better for GVA2 if m > ni, as will typically be the case in 
  applications. When m < ni, GVA seems to perform better. We're still faster in both cases
  though.
- Look at Journal of Computational Statistics
- Work on presentation on tools:
  * GitHub/Version control
  * Package development
  * Profiling and debugging
  * AWS
  * Text editors/Sublime Text as an example
  * UNIX
  * Hadley Wickham's various contributions to the R ecosystem: testthat, magrittr, dplyr
  * Twitter and the blogosphere
  * Leechblock
- Application to physical activity data
- Read more about numerical linear algebra and optimisation
- Think about speeding things up with mclapply
  Done. This was really easy. I should do the same thing with coverage_percentage.
- Fill in Tutorial Preferences form
- Drop off solutions to Diana Warren
  Done
- Change convergence test in zero_inflated_model to check whether lower bound has stopped
  increases _and_ parameters have stopped changing.
- What am I going to write about in my paper?
  Look at draft
  Look up expressions for accuracy of solving systems, so that you can compare methods of
  solving \mLambda and \mLambda^-1

15/06/2015
----------
- Error analysis for solving linear systems
  https://www.fields.utoronto.ca/programs/scientific/09-10/FoCM/courses/Cucker5.pdf
- Rewrote median_accuracy.R as a script that accepts command line arguments, so that you
  can run it from the command line.
- Found another overflow in Laplace's method in gaussian.R

04/06/2015
----------
- Tutored last Thurs 10am class
  They were a little more attentive, and I didn't make any arithmetic mistakes.
- Attended Winter School MATH1005 meeting
- Reading papers on zero-inflated models

03/06/2015
----------
- Tutored last Wed 11am class
  Only half the class was there. Most of them are already in holiday mode. Some people asked
  good exam questions.
- Spent more time trying alternate expressions for the derivative of GVA2.
- Got chewed out by Samuel for not submitting my Postgraduate School Progress Report and
  Training Needs Analysis
  I've completed it now, and just need John to sign off on his part.
- Run simulations
  I'm running the 1e6 samples MCMCs for random intercepts and random slopes.
- Median accuracy/Coverage?
- Literature review
- Read more numerical linear algebra/numerical analysis
  James Demmel's books and lecture notes are really good
- Prepare presentation

31/05/2015
----------
- Install SublimeREPL and R-Box packages in Sublime. Then you can work in Sublime most of the time.
  Done.
- Install Monokai-Black, just so you feel at home.
- Tried Matrix library again.
  mR <- new("dtrMatrix", uplo="L", diag="N", x=as.vector(decode$mR), Dim=as.integer(c(d, d)))
  Just as disappointing as always. Constructing the matrix takes more time than using it 
  saves.
- I got fed up with pushing against the R performance wall, and wrote fastsolve in C++ using
  Eigen.

28/05/2015
----------
- Fixed overflow in GVA (!)
- Derivatives in GVA and GVA2 are correct, checked numerically
- All accuracy tests run

25/05/2015
----------
- Finish draft of Progress Report, and send it to Samuel
- Continue investigating matrix derivative in GVA algorithm
  The matrix derivative seems definitely wrong.
  I have some ideas as to why, but my attempts to fix it have just left me more confused.
- The use of vtheta in various places in the code base is getting in my way. I'm
  going to make sure that vtheta is only a parameter in functions which get called directly by
  optim.
- Prepare presentation of progress
  Show spline fits
  Show effective sample sizes

  The meeting was cancelled today

21/05/2015
----------
- GVA succeeds on the spline test case after five iterations
- GVA2 currently fails with a non-finite value supplied by optim

18/05/2015
----------
- John said that BSpline models can be non-identifiable.
  * I should check the chains, and
  * plot the functions based on the coefficients that I've fit.
  How do I plot the functions? I have a vector of co-efficients, but haven't yet found out
  how to evaluate the BSpline functions using them
- Finish the report for Samuel Mueller.
- Fill in roles.
- Report missing mark for student.

14/05/2015
----------
- Returned the ute
- Tutored
- Begun filling in Postgraduate School Progress and Training Needs Analysis
  It's complicated - I need to summarise my progress to date and training needs. So I'll
  write a draft first, and then write my final copy.

13/05/2015
----------
- Tutored
- Hired a ute to transport Dex's stuff
- Re-assembled Dex's loft bed

12//05/2015
-----------
- Disassembled Dex's loft bed

11/05/2015
----------
- I prepared a presentation to update people on my progress. I asked Shila for some
  time, which she allocated. But in the end, we ran out of time.

  I sent the presentation and associated graphs to John.
- I'm supposed to present in two weeks. I don't expect to be able to complete anything
  major in that time, but hopefully I can make enough progress that I have something to
  show.
- Susan Lin said that she got two out of three wrong in this week's quiz, and she wondered
  if there was a mistake.
- Pay claims, assignment marks
- Training Report for Samuel

07/05/2015
----------
- Write up meeting notes
  Done
- Use John's new formOmega code to construct mZ
- Fit model
  Almost worked!
- Fill in report for Samuel

05/05/2015
----------
- Get random slopes working again
  Done. It works quite nicely now.
- Get splines working.
  Fail. Nothing I'm trying works. I don't know what the problem is.
- Correct derivative in GVA code
  I still don't believe this is right. Check with numeric differentiation.
- Attend Advanced R reading group
  Done

04/05/2015
----------
- Attend meeting
  We talked about model selection
- Continue to debug spline problems
  * ZOSull doesn't seem to work properly in the range [-pi/2, pi/2]. Works much better in
    the range [10, 100].
  * Stan always seems to get roughly the same answer, even if I change the function, whereas
    the variational approximation gives different answers. This makes me suspect that Stan
    is doing something strange.

30/04/2015
----------
- Give tutorial
  Done
- I worked out that the function I was trying to fit looked like noise. I changed it to
  something more reasonable.
- Caught a bug in the mean field update. I was using block size in the divisor for the Inverse
  Wishart update which is correct when doing random intercepts or random slopes, but wrong when 
  doing splines

29/04/2015
----------
- Prepare for tutorial
  Done, but I was cutting it fine
- Give tutorial
  Done
- Mark assignments - due next week
- Continue modifying the software to work p=0
  It's producing garbage results at present
- Fix bugs in the lower bound
- Try to rationalise the notation in the draft of the paper I've done
  There were two scalars both named p, one the dimension of the fixed effects and the other
  the dimension of the inverse Wishart distribution that we're using. This was quite a mess.
  I've introduced a new variable d to play the wrole of the second scalar.
- Something funny is going on in ZOSull once we calculate Z = B %*% LZ
  B is nice and banded as we expect, but Z = B LZ is mostly small values/zeroes.


28/04/2015
----------
- Sidestep repeated columns problem by building a spline model that uses random effects
  only
  I started modifying the software to be allow no fixed effects in a model. This turned up
  all sorts of errors and edge cases in the code which needed to be dispensed with.
- Attend Advanced R reading group

27/04/2015
----------
- Spent an hour and a half understanding ten lines of code that calculates the Omega
  penalty matrix in ZOSull.

- Install R packages
To install limma:
source("http://bioconductor.org/biocLite.R")
biocLite("limma")

To install rstan:
source('http://mc-stan.org/rstan/install.R', echo = TRUE, max.deparse.length = 2000)
install_rstan()

To install the other packages from CRAN:
install.packages(c("Rcpp", "RcppEigen"))

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

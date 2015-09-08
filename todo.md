To Do list
==========

26/08/2015
----------
- Continue writing up application
  Most numeric results finished. Now need accuracy.
- Meet with John
- Computer practical

25/08/2015
----------
- Write about application
- Give draft to John
- Fix problem in application
  Done - it was because the senior column in mX was highly collinear with mZ. I have no idea
  why this was true, but I don't even really care about the senior column. I was just including
  it to adjust for more covariates.
- Accuracy for application
- What graphs do you need?
- Think about getting an nVidia card
- Moving
- Vaccinations

24/08/2015
----------
- Passport
  Sent! Wow, I actually did it.
- Prepare for teaching
  Done
- Continue with application
  It runs, but there's some weird stuff in diag(mLambda).
- Bayesian model selection
- Write!
  You could draft the section on the dimension reduction
- Claim pay
  Done for 27/7, 28/7, 29/7, 3/8, 4/8, 5/8, 10/8, 11/8, 12/8, 17/8, 18/8, 19/8
  Not sure what to do about Winter School. I'm afraid that Sonia's going to scream at me.
- Pack

19/08/2015
----------
- Try new dimension reduction code on the application example.
  It worked, but from my point of view, the results are disappointing.
  Time difference of 2.963303 mins
  Time difference of 55.16225 secs  
- Teach
  Done
- Enter marks for last week's computer practicals
  Done
- Claim money
- Find data set
- Pack
- Passport

18/08/2015
----------
- I tried the dimension reduction idea. So far, it's not panning out. But I wonder if I need
  to re-order the covariance matrix. I still have code for that in the repository.

  It worked!

14/08/2015
----------
- Why doesn't coverage percentage work anymore?
  Why is E[\vbeta|\vy, \mX] not what you think it is?
  log(\vr \vy) =  (\mX \vbeta + \mZ \vu)
  Is it the test case?
  Could you get the expected parameters from the MCMC runs?
  If I could just get that to work ...
- PH hasn't responded to my data request, and she might not.
  Need a backup plan.
- Speed
  Need to nail this down once and for all. I think that the L-BFGS-B traversal is faster for
  GVA2 than for GVA. But I've made a lot of changes since I tested this last, which makes me nervous.
  Still seems to be true. Running accuracy.R -t slope I get:
  Time difference of 2.472427 secs
  [1] 0.6580959 0.6919876
  [1] 0.7191744
  [1] 0.5296653 0.5219726
  [1] 0.9107536
  Time difference of 4.899895 secs
  [1] 0.8786032 0.8833411
  [1] 0.9080758
  [1] 0.5235579 0.5155404
  [1] 0.8958491
  Time difference of 2.40483 secs
  [1] 0.8932545 0.8988409
  [1] 0.9117547
  [1] 0.5235019 0.5154443
  [1] 0.8958394
  Time difference of 0.4393022 secs
  [1] 0.8786032 0.8833409
  [1] 0.9080757
  [1] 0.5235580 0.5155404
  [1] 0.8958491

  markg@markg-OptiPlex-9020:~/phd/code$ grep -i iteration /tmp/gva_log.txt | wc -l
  3810
  markg@markg-OptiPlex-9020:~/phd/code$ grep -i iteration /tmp/gva2_log.txt | wc -l
  849  
  
  That's that.
  
- John wants me to take advantage of the sparsity of \mLambda = (\mR \mR)^\top by making
  Rinds include only the nonzero entries of \mLambda. I really doubt this is going to help very
  much, as it will only effect the execution of optim(), but we'll see. It's worth a try.
  mC2 <- mC
  mC2[mC != 0] <- 1
  ones <- crossprod(mC2) + solve(mSigma.inv)
  Rinds <- lower.tri(mLambda) intersect ones

13/08/2015
----------
- Why doesn't coverage percentage work anymore?
  Why is E[\vbeta|\vy, \mX] not what you think it is?
- Fix up tickets with Jetstar
  Done
- Passport application
  Entered details, form printed, paid for.
- Chase Ben about boxes
- Book car for Lewys
- Lecture
  Done
- Pay claim
  Winter School - What days did I work? That's quite a lot of money.
  Tutorial so far this semester. That's also quite a lot of money.
- Median graph
  Done, although not very pretty
- 

12/08/2015
----------
- Email PH about data release
  Done
- Email HR with passport information
  Done
- Enter marks for computer prac last week
  Done
- Renew British passport
- Give computer prac at 3pm
- Entered some numerical results for variance components of random slopes model in
  paper draft
- Coverage percentiles
  Running
- Figure out what to do about graphs you need to produce, and where in your draft paper
  they need to go
  Need to produce median accuracy graphs which have different kinds of approximations side by
  side, with correct labelling

10/08/2015
----------
- Prepare for STAT2012 tutorials
  Done. That didn't go very well, I'll prepare better for tomorrow.
- I was able to get the application to work last week, so it's time to move on.
  Gather more results
  Graphs
  Write
- Replace B() functions in draft paper with explicit functions
- Email PRC about using the data
- Think about thesis
  Numeric problems with GVA
  Optimisation/domain of attraction
  Zero inflated Poisson models
- Mark lab reports from STAT2012
- Claim pay for Winter School
- Claim pay for work to date during this semester
- Send scan of passport to HR
- Fix sigma accuracy calculation.
  I'm really not looking forward to this.
- Coverage percentile

03/08/2015
----------
- Application
  I'm looking at the PRC data set. There's a sea of variables.
    [1] "Id"                      "age"                     "agegrp"                 
    [4] "Gender"                  "Education"               "Employment"             
    [7] "Income"                  "Language"                "Indigenous"             
   [10] "enough_Income"           "SEIFA_quintile"          "aria_mean"              
   [13] "aria_cat"                "Gprequired"              "suff_pa"                
   [16] "enough_fruit"            "enough_veg"              "waistcircrisk"          
   [19] "baselinedata"            "threemonth"              "twelvemonth"            
   [22] "heightcm"                "time"                    "weight"                 
   [25] "waistcirc"               "fruitday"                "vegday"                 
   [28] "chips"                   "redmeat"                 "sausage"                
   [31] "takeaway"                "nuts"                    "cakes"                  
   [34] "sweets"                  "legumes"                 "cereal"                 
   [37] "trimfatmeat"             "olive"                   "otheroil"               
   [40] "margarine"               "butter"                  "chickentrim"            
   [43] "lowfmilk"                "litecream"               "lfcheese"               
   [46] "spreadlowf"              "wholepasta"              "brownrice"              
   [49] "breads"                  "PAtired"                 "PAstress"               
   [52] "PAdemands"               "PAdepress"               "diethurry"              
   [55] "dietstress"              "dietdemands"             "dietdepress"            
   [58] "oftenPA"                 "PAwhen"                  "wherePA"                
   [61] "whomPA"                  "rainPA"                  "goalPA"                 
   [64] "incidgoal"               "whencycle"               "oftenwalk"              
   [67] "wupstairs"               "incidbackup"             "cookbook"               
   [70] "functiona"               "tablenice"               "planhealthy"            
   [73] "preparef"                "fameasier"               "PAwithyou"              
   [76] "regPAfam"                "choicewgoal"             "PAhelpgoal"             
   [79] "sadcheer"                "nervous"                 "restless"               
   [82] "hopeless"                "effort"                  "worthless"              
   [85] "timeswalked"             "modtimes"                "vigtimes"               
   [88] "minwalked"               "minvig"                  "bmi_cat"                
   [91] "bmi"                     "too_short"               "tottimes"               
   [94] "mintot"                  "minmod"                  "fat_score"              
   [97] "fibre_score"             "total_score"             "social_support_pa"      
  [100] "social_support_eating"   "planning_healthy_eating" "planning_incidental_pa" 
  [103] "planning_structured_pa"  "confidence_eating"       "confidence_pa"

  Splines, mixed models, what to do?
  Can you find and read the paper? Or just remember what was important?
  What should go into mX and mZ? What should the covariance matrix look like?

  Update: Okay, I've thought about it a bit. The primary outcome of interest is the positive
  relationship between time and total minutes of physical activity aka mintot. There are
  273 individuals, observed at 3 time points, 0, 3 and 12 months.
  A random slope model on time with zero inflation seems appropriate, using fixed effects to
  adjust for other covariates of interest. Adjust for all the standard things, socioeconomic
  and demographic factors, and of course age.
   
30/07/2015
----------
- Spline accuracy results are now good (!)
- Application
- Marking

29/07/2015
----------
- Try fixing spline accuracy problems by lowering the maximum number of VB iterations
  to 20.
  The problem turned out to be the MCMC! Generating a different data set from the same
  seed and re-running fixed it.

28/07/2015
----------
- Poor accuracy for Test case 97 has something to do with diag(mLambda) having some large
  values.
  > vaccuracy
    [1] 0.8970096 0.8539605 0.9023597 0.8229834 0.9060690 0.8744393 0.8984698
    [8] 0.8354486 0.8441881 0.8366873 0.7080212 0.8499089 0.9160671 0.8074408
   [15] 0.8082496 0.7628658 0.8706933 0.8767969 0.7005117 0.8522690 0.8163083
   [22] 0.8784463 0.7271954 0.6014027 0.4452374 0.7133955 0.8697347 0.8842829
   [29] 0.8028537 0.7025526 0.5527570 0.9047466 0.8384766 0.8831617 0.8355878
   [36] 0.8853087 0.6911480 0.8992032 0.8451538 0.9117479 0.8412023 0.8538672
   [43] 0.8649741 0.8747038 0.5066782 0.8723054 0.8013161 0.8542742 0.8870647
   [50] 0.7792468 0.8829639 0.4386193 0.8283087 0.7846913 0.8583867 0.8405602
   [57] 0.8241351 0.7029739 0.9174753 0.5346066 0.5671533 0.7324143 0.8498585
   [64] 0.5499912 0.7246882 0.8484476 0.8996635 0.8119486 0.7715703 0.9362681
   [71] 0.6509210 0.7006210 0.6682824 0.8349544 0.7625612 0.6258402 0.8863403
   [78] 0.6876932 0.6964111 0.8551330 0.8045416 0.8036846 0.7262400 0.7989986
   [85] 0.6739237 0.8345746 0.6920244 0.8490511 0.8819539 0.6701659 0.9095995
   [92] 0.8956516 0.7851022 0.8487144 0.7899954 0.8580233 0.4214853 0.9157717
   [99] 0.8129945 0.9029803
  > order(vaccuracy)
    [1]  97  52  25  45  60  64  31  61  24  76  71  73  90  85  78  37  87  79
   [19]  19  72  30  58  11  26  65  83  23  62  75  16  69  50  54  93  95  84
   [37]  47  29  82  81  14  15  68  99  21   4  57  53  86  74   8  35  10  33
   [55]  56  41   9  39  66  94  88  63  12  20  42   2  48  80  96  55  43  27
   [73]  17  46   6  44  18  22  89  51  34  28  36  77  49  92   1   7  38  67
   [91]   3 100  32   5  91  40  98  13  59  70  

   Having looked at the logs, I think the problem is that if we keep iterating
   beyond a certain point we go from a good solution to a ludicrous one. This is
   because as vmu goes to silly negative values, the diagonals of mLambda become
   large to accomodate this. I'm trying to side step the problem by doing less
   iterations.

   My computation is taking forever to run today.
- Prepare for two o'clock tutorial.
- Application

27/07/2015
----------
- Check accuracy results for splines for GVA
  Done. There are bad results for a small number of test cases: 
- Run median accuracy for splines for other approximations
  Done.
- Pay claims for Winter School
- Prepare for tutoring second year stats
  Jennifer did it instead.
- Prepare presentation for melanoma meeting
  Done
- Invoice for Meanpath
- Proof of citizenship for University of Sydney
- Submit results from Winter School to Dianne
  Done

22/07/2015
----------
- Spent the last few days converting one of John's R programs into C++. This resulted in
  a 13.3 times speed-up (!).
- Package the C++ code up in R
- Spline median accuracy, focus on the functions themselves
  It should just be accuracy(true function, fitted function) across all MCMC runs and
  approximations
  i.e. loop through all the data sets, loading them
  do the VB approximation
  note the accuracy
  graph it
- Claim pay
  Last semester
  Winter school

15/06/2015
----------
- GVA2
  real  22m32.521s
  user  22m19.472s
  sys 0m8.429s

  Will GVA be faster? Or slower?
  real  45m52.144s
  user  45m22.574s
  sys 0m12.253s

  Quite a difference.

- Splines next
  This is going to be complicated, I expect everything will break.

08/06/2015
----------
- Re-run median accuracy for slope test case with the new MCMC data
  GVA crashed on test case 4. This turned out to be due to normalisation issues in the
  accuracy integral. So I've fixed that.
  GVA and GVA2 work pretty reliably. GVA2 is a little faster.
  GVA_nr is very fast, but unstable/unreliable. Occasionally the lower bound becomes infinite
  and the results are complete nonsense.
  Laplace takes a very long time to converge, and the accuracy results aren't that great.

07/06/2015
----------
- Re-run MCMC
  I'd hardcoded some settings. That came back to bite me. I have to re-run the whole
  task now.
- R linting
- Writing
- Check quiz answers

06/06/2015
----------
- Fix crashes coming from integrate() in median_accuracy()
- Fix crash coming from create_accuracy_df() in median_accuracy_graph()
  Done. We can now produce a median accuracy graph.
- Mark reports
  Done
- Fill in spreadsheet of reports
  Done
- Fill in HR
- Meanpath invoice
- Scan report marking template sheets
- Respond to Suzanne's email re: True stories

02/06/2015
----------
- In random slopes median accuracy, occasionally mLambda's diagonal contains a huge
  value and the rest of the fit is garbage. T1 in the lower bound is usually Inf in these
  cases.

30/06/2015
----------
- Found major bug:
    mSigma.u.inv <- kronecker(diag(1, m - 1), solve(mPsi)/(v - d - 1))
  should have been
    mSigma.u.inv <- kronecker(diag(1, m - 1), solve(mPsi/(v - d - 1)))
  That is, I was inverting and then dividing, rather than dividing and then inverting.

  This had the effect of raising the accuracy of _all_ the intercepts above 80%!
- To do: mean of random effects
  Re-arrange boxplots so that accuracy of methods can be compared
  Think of revisions to paper
  Median accuracy for slopes and splines
  
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

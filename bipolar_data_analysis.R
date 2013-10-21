# bipolar_data_analysis.R
setwd("~/Desktop/cchs2002/DATA")
hs = read.csv("hs.csv.gz")
dim(hs)
require(Hmisc)
sink("~/phd/ccsh_describe.txt")
describe(hs)
sink()
# Make names lowercase, because I find this much easier to type and to
# read
names(hs) = tolower(names(hs))
find_column = function(pattern) {
  return(grep(pattern, names(hs), ignore.case = TRUE, value = TRUE))
}
find_column("age")
find_column("dep.*28")

attach(hs)
# depressive = dep_q38?
table(depb_38)
# mania
table(miab_29)
# Could classify anything over 0 as a yes. But what do we do about the don't
# knows, the not applicables and the not stateds?

# Table 1 - sociodemographics ----
table(dhhb_sex)
# Q: Do we have the continuous age variable? From the analysis in Table 1,
# our collaborators must have this.
# Age
table(dhhbgage)
# Education
table(edubdr04)
# Employment
table(lbfbg31)

# Table 1 - clinical ----
# Suicide attempts 12 month
table(depbfsya)
# Number of days unable to work
table(depb_68)
# Interf. with work/social/personal
# Depression/mania interfered with work/school/personal?
# miabdint/padq_22 for mania? Which variable did they use?
# Depressive feelings interfere with work/social/relationship - depb_28
table(depb_28)
# Which variable did they use for Medication 12 mo? There are several
# variables that look relevant, but I'm not sure which one to use.

detach()
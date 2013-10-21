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

# Table 2 Summary of social support by mood disorder ----
# Mean number of friends and relatives
summary(as.numeric(ssmb_01))
# Summary of types of social support
# Tangible social support
table(ssmb_21a)
# Affective support
table(ssmb_22a)
# Positive social interactive support
table(ssmb_23a)
# Emotional/informational support
table(ssmb_24a)

# Table 3 Prevalence of specific support by mood disorders ----
# Tangible social support
# Has someone to give help if confined to bed
table(ssmb_02)
# Has someone to take to doctor
table(ssmb_05)
# Has someone to prepare meals
table(ssmb_12)
# Has someone to help with daily chores if sick
table(ssmb_15)

# Emotional or informational support
# Has someone to listen
table(ssmb_03)
# Has someone to give advice about a crisis
table(ssmb_04)
# Has someone who gives information to understand a situation
table(ssmb_08)
# Has someone to confide in
table(ssmb_09)
# Has someone to share most private worries and fears with
table(ssmb_16)
# Has someone to turn to for suggestions for personal problems
table(ssmb_17)
# Has someone who understands problems
table(ssmb_19)

# Positive social interaction
# Has someone to have a good time with
table(ssmb_07)
# Has someone to get together with for relaxation
table(ssmb_11)
# Has someone to get mind off things
table(ssmb_14)
# Has someone to do something enjoyable with
table(ssmb_18)

# Affective support
# Has someone who shows love and affection
table(ssmb_06)
# Has somsone who gives hugs
table(ssmb_10)
# Has someone who loves and makes feel wanted
table(ssmb_20)

detach()


# bipolar_data_analysis.R
require(Hmisc)
require(survey)
require(foreign)
#require(sqldf)

# Load the data, and create a survey object ----
#setwd("~/Documents/cchs2002/DATA")
hs = read.csv("hs.csv.gz")
names(hs) = tolower(names(hs))
cchs = read.spss("Cchs1-2.sav coping Aug 25 2013.sav")
cchs = as.data.frame(cchs)
dim(hs)
dim(cchs)
intersect(names(hs), names(cchs))
setdiff(names(cchs), names(hs))
names(hs) = tolower(names(hs))
#hs_joined = sqldf("select a.*, b.Mood12, b.Bip12, b.Dep12 from hs a left join cchs b on a.admb_rno=b.admb_rno")
hs_joined = merge(hs, cchs[,c("admb_rno", "Mood12", "Bip12", "Dep12")], by = "admb_rno")
hssvy = svydesign(~1, weights=~wtsb_m, data=hs_joined)

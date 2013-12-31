# bipolar_data_analysis.R
require(Hmisc)
require(survey)
require(foreign)
#require(sqldf)

# Load the data, and create a survey object ----
setwd("~/Documents/cchs2002/DATA")
hs = read.csv("hs.csv.gz")
names(hs) = tolower(names(hs))
# Clean and prepare data ----
for (i in 1:nrow(hs)) {
  row = hs[i,]
  
  if (row$miabdey == 1)
    row$deponly = 0
  else if (row$depbddy==1)
    row$deponly = 1
  else
    row$deponly = NA
  
  if (row$miabdey == 1)
    row$mood = 0
  else if (row$depbddy==1)
    row$mood = 1
  else
    row$mood = NA
  
  if (row$aldb_01==1 | row$aldb_03==1 | (7 <= row$aldb_15c && row$aldb_15c <= 10) | (7 <= row$aldb_15d && row$aldb_15d<= 10)) {
    row$alc_abuse = 1
  }
  else if (row$aldb_01==2 & row$aldb_03==2 & (0 <= row$aldb_15c & row$aldb_15c <= 6) & (0 <= row$aldb_15d & row$aldb_15d <= 6)) {    
    row$alc_abuse = 0
  } else {
    row$alc_abuse = NA
  }
  
  set_missing = function(row, fields, threshold_type, threshold)
  {
    stopifnot(threshold_type %in% c("equals", "greater-than"))
    for (field in fields) {
      if (threshold_type == "equals") {
        if (row[field] == threshold) {
          row[field] = NA
        }
      }
      if (threshold_type == "greater-than") {
        if (row[field] >= threshold) {
          row[field] = NA
        }
      }
      
    }
    return(row)
  }

  row = set_missing(row, c("depbddy", "depbfsya", "edubdh04", "depbfsla", "padbddy", "sopbdpy"), "equals", 9)
  row = set_missing(row, c("incb_4"), "greater-than", 999996)  
  row = set_missing(row, c("lbfb_01", "depb_28", "medb_11c",  "medb_11d", "medb_11e",
             "medb_11f", "miab_09", "ssmb_21a", "ssmb_22a", "ssmb_23a", "ssmb_24a", "aldb_01", "aldb_03"), "greater-than", 6)

  row = set_missing(row, c("depbdon", "miabdon", "depb_68", "depb_721", "depb_871", "miab_29", "miab_33a", "miab_48a"), "greater-than", 996)
  row = set_missing(row, c("depb_372", "miab_18c"), "greater-than", 3)
  row = set_missing(row, c("depb_88", "miab_481", "aldb_15c", "aldb_15d"), "greater-than", 96)
  row = set_missing(row, c("aldbint"), "greater-than", 99.6)
  row = set_missing(row, c("aldbdpp"), "greater-than", 9.99)
  row = set_missing(row, c("ssmb_01"), "greater-than", 997)
  row = set_missing(row, c("ssmb_02", "ssmb_03", "ssmb_04", "ssmb_05", "ssmb_06", "ssmb_07", "ssmb_08", "ssmb_09", "ssmb_10",
                           "ssmb_11", "ssmb_12", "ssmb_13", "ssmb_14", "ssmb_15", "ssmb_16", "ssmb_17", "ssmb_18", "ssmb_19", "ssmb_20",
                           "strb_65b", "strb_65c", "strb_65d", "strb_66", "strb_67"), "greater-than", 7)
  if (row$strb_65b == 5)
    row$strb_65b = NA
  
  if (row$aldbdsf == 99) {
    row$aldbdsf = NA
  }
  
  # Create factors according to email from Dr Wang ----
  if (row$strb_63 == 1 | row$strb_64 == 1 | row$strb_65a == 1 | row$strb_610 == 1 | row$strb_611 == 1)
    row$avoidance = 1
  else if (any(is.na(row[c("strb_63", "strb_64", "strb_65a", "strb_610", "strb_611")])))
    row$avoidance = NA
  else
    row$avoidance = 0
  
  if (row$strb_61 %in% c(3, 4) | row$strb_62 %in% c(3, 4) | row$strb_68 %in% c(3, 4) | row$strb_69 %in% c(3, 4))
    row$prob_solving = 1
  else if (any(is.na(row[c("strb_61", "strb_62", "strb_68", "strb_69")])))
    row$prob_solving = NA
  else
    row$prob_solving = 0
  
  if (row$strb_65b == 1 | row$strb_65c == 1 | row$strb_65d == 1 | row$strb_66 %in% c(3, 4) |
        row$strb_67 %in% c(3, 4)) {
    row$behaviours = 1
  }
  else if (any(is.na(row[c("strb_65b", "strb_65c", "strb_65d", "strb_66", "strb_67")]))) {
    row$behaviours = NA
  }
  else
    row$behaviours = 0
}
mood = subset(hs, !is.na(deponly))
hs = mood

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

find_column = function(pattern) {
  return(grep(pattern, names(hs), ignore.case = TRUE, value = TRUE))
}

# New analysis - started from scratch ----


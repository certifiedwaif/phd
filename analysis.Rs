####################################################################################################

#source("mergeData.Rs")
source("stepwise.r")

all.na <- function(x) { return( all(is.na(x)) ); }
any.na <- function(x) { return( any(is.na(x)) ); }
sum.na.rm <- function(x) { sum(x,na.rm=TRUE); }

MIN_OBS_IN_CATEGORY <- 20

####################################################################################################

# Select subset of the data
#cchs_mod  <- cchs[cchs$Mood12!="control",]
cchs_mod<- cchs[cchs$Mood12=="bipolar",]
#cchs_mod<- cchs[cchs$Mood12=="depression",]
#cchs_mod<- cchs 

####################################################################################################
	
# Remove any rows that contain all missing values
cchs_mod <- cchs_mod[!apply(cchs_mod,1,all.na),]

# Remove any columns that contain all missing values
cchs_mod <- cchs_mod[,!apply(cchs_mod,2,all.na)]

####################################################################################################

# Creat some derived variables

n <- length(cchs_mod$strb_61)
	
strb_61_num <- rep(NA,n)
strb_61_num[which(cchs_mod$strb_61=="OFTEN")] <- 3
strb_61_num[which(cchs_mod$strb_61=="SOMETIMES")] <- 2
strb_61_num[which(cchs_mod$strb_61=="RARELY")] <- 1
strb_61_num[which(cchs_mod$strb_61=="NEVER")] <- 0

strb_61_OFTEN 	  <- (cchs_mod$strb_61=="OFTEN")
strb_61_SOMETIMES <- (cchs_mod$strb_61=="SOMETIMES")
strb_61_RARELY    <- (cchs_mod$strb_61=="RARELY")
strb_61_YES       <- (cchs_mod$strb_61=="OFTEN")|(cchs_mod$strb_61=="SOMETIMES")
	
strb_62_num <- rep(NA,n)
strb_62_num[which(cchs_mod$strb_62=="OFTEN")] <- 3
strb_62_num[which(cchs_mod$strb_62=="SOMETIMES")] <- 2
strb_62_num[which(cchs_mod$strb_62=="RARELY")] <- 1
strb_62_num[which(cchs_mod$strb_62=="NEVER")] <- 0

strb_62_OFTEN 	  <- (cchs_mod$strb_62=="OFTEN")
strb_62_SOMETIMES <- (cchs_mod$strb_62=="SOMETIMES")
strb_62_RARELY    <- (cchs_mod$strb_62=="RARELY")
strb_62_YES       <- (cchs_mod$strb_62=="OFTEN")|(cchs_mod$strb_62=="SOMETIMES")
	
strb_63_num <- rep(NA,n)
strb_63_num[which(cchs_mod$strb_63=="OFTEN")] <- 3
strb_63_num[which(cchs_mod$strb_63=="SOMETIMES")] <- 2
strb_63_num[which(cchs_mod$strb_63=="RARELY")] <- 1
strb_63_num[which(cchs_mod$strb_63=="NEVER")] <- 0

strb_63_OFTEN 	  <- (cchs_mod$strb_63=="OFTEN")
strb_63_SOMETIMES <- (cchs_mod$strb_63=="SOMETIMES")
strb_63_RARELY    <- (cchs_mod$strb_63=="RARELY")
strb_63_YES       <- (cchs_mod$strb_63=="OFTEN")|(cchs_mod$strb_63=="SOMETIMES")
	
strb_64_num <- rep(NA,n)
strb_64_num[which(cchs_mod$strb_64=="OFTEN")] <- 3
strb_64_num[which(cchs_mod$strb_64=="SOMETIMES")] <- 2
strb_64_num[which(cchs_mod$strb_64=="RARELY")] <- 1
strb_64_num[which(cchs_mod$strb_64=="NEVER")] <- 0

strb_64_OFTEN 	  <- (cchs_mod$strb_64=="OFTEN")
strb_64_SOMETIMES <- (cchs_mod$strb_64=="SOMETIMES")
strb_64_RARELY    <- (cchs_mod$strb_64=="RARELY")
strb_64_YES       <- (cchs_mod$strb_64=="OFTEN")|(cchs_mod$strb_64=="SOMETIMES")
	
strb_65A_num <- rep(NA,n)
strb_65A_num[which(cchs_mod$strb_65a=="OFTEN")] <- 3
strb_65A_num[which(cchs_mod$strb_65a=="SOMETIMES")] <- 2
strb_65A_num[which(cchs_mod$strb_65a=="RARELY")] <- 1
strb_65A_num[which(cchs_mod$strb_65a=="NEVER")] <- 0

strb_65A_OFTEN 	   <- (cchs_mod$strb_65a=="OFTEN")
strb_65A_SOMETIMES <- (cchs_mod$strb_65a=="SOMETIMES")
strb_65A_RARELY    <- (cchs_mod$strb_65a=="RARELY")
strb_65A_YES       <- (cchs_mod$strb_65a=="OFTEN")|(cchs_mod$strb_65a=="SOMETIMES")
	
strb_65B_num <- rep(NA,n)
strb_65B_num[which(cchs_mod$strb_65b=="OFTEN")] <- 3
strb_65B_num[which(cchs_mod$strb_65b=="SOMETIMES")] <- 2
strb_65B_num[which(cchs_mod$strb_65b=="RARELY")] <- 1
strb_65B_num[which(cchs_mod$strb_65b=="NEVER")] <- 0
strb_65B_num[which(cchs_mod$strb_65b=="DO NOT SMOKE")] <- 0

strb_65B_OFTEN 	   <- (cchs_mod$strb_65b=="OFTEN")
strb_65B_SOMETIMES <- (cchs_mod$strb_65b=="SOMETIMES")
strb_65B_RARELY    <- (cchs_mod$strb_65b=="RARELY")
strb_65B_YES       <- (cchs_mod$strb_65b=="OFTEN")|(cchs_mod$strb_65b =="SOMETIMES")
	
strb_65C_num <- rep(NA,n)
strb_65C_num[which(cchs_mod$strb_65c=="OFTEN")] <- 3
strb_65C_num[which(cchs_mod$strb_65c=="SOMETIMES")] <- 2
strb_65C_num[which(cchs_mod$strb_65c=="RARELY")] <- 1
strb_65C_num[which(cchs_mod$strb_65c=="NEVER")] <- 0

strb_65C_OFTEN 	   <- (cchs_mod$strb_65c=="OFTEN")
strb_65C_SOMETIMES <- (cchs_mod$strb_65c=="SOMETIMES")
strb_65C_RARELY    <- (cchs_mod$strb_65c=="RARELY")
strb_65C_YES       <- (cchs_mod$strb_65c=="OFTEN")|(cchs_mod$strb_65c =="SOMETIMES")
	
strb_65D_num <- rep(NA,n)
strb_65D_num[which(cchs_mod$strb_65d=="OFTEN")] <- 3
strb_65D_num[which(cchs_mod$strb_65d=="SOMETIMES")] <- 2
strb_65D_num[which(cchs_mod$strb_65d=="RARELY")] <- 1
strb_65D_num[which(cchs_mod$strb_65d=="NEVER")] <- 0

strb_65D_OFTEN 	   <- (cchs_mod$strb_65d=="OFTEN")
strb_65D_SOMETIMES <- (cchs_mod$strb_65d=="SOMETIMES")
strb_65D_RARELY    <- (cchs_mod$strb_65d=="RARELY")
strb_65D_YES       <- (cchs_mod$strb_65d=="OFTEN")|(cchs_mod$strb_65d =="SOMETIMES")
	
strb_66_num <- rep(NA,n)
strb_66_num[which(cchs_mod$strb_66=="OFTEN")] <- 3
strb_66_num[which(cchs_mod$strb_66=="SOMETIMES")] <- 2
strb_66_num[which(cchs_mod$strb_66=="RARELY")] <- 1
strb_66_num[which(cchs_mod$strb_66=="NEVER")] <- 0

strb_66_OFTEN 	  <- (cchs_mod$strb_66=="OFTEN")
strb_66_SOMETIMES <- (cchs_mod$strb_66=="SOMETIMES")
strb_66_RARELY    <- (cchs_mod$strb_66=="RARELY")
strb_66_YES       <- (cchs_mod$strb_66=="OFTEN")|(cchs_mod$strb_66 =="SOMETIMES")
	
strb_67_num <- rep(NA,n)
strb_67_num[which(cchs_mod$strb_67=="OFTEN")] <- 3
strb_67_num[which(cchs_mod$strb_67=="SOMETIMES")] <- 2
strb_67_num[which(cchs_mod$strb_67=="RARELY")] <- 1
strb_67_num[which(cchs_mod$strb_67=="NEVER")] <- 0

strb_67_OFTEN 	  <- (cchs_mod$strb_67=="OFTEN")
strb_67_SOMETIMES <- (cchs_mod$strb_67=="SOMETIMES")
strb_67_RARELY    <- (cchs_mod$strb_67=="RARELY")
strb_67_YES       <- (cchs_mod$strb_67=="OFTEN")|(cchs_mod$strb_67=="SOMETIMES")
	
strb_68_num <- rep(NA,n)
strb_68_num[which(cchs_mod$strb_68=="OFTEN")] <- 3
strb_68_num[which(cchs_mod$strb_68=="SOMETIMES")] <- 2
strb_68_num[which(cchs_mod$strb_68=="RARELY")] <- 1
strb_68_num[which(cchs_mod$strb_68=="NEVER")] <- 0

strb_68_OFTEN 	  <- (cchs_mod$strb_68=="OFTEN")
strb_68_SOMETIMES <- (cchs_mod$strb_68=="SOMETIMES")
strb_68_RARELY    <- (cchs_mod$strb_68=="RARELY")
strb_68_YES       <- (cchs_mod$strb_68=="OFTEN")|(cchs_mod$strb_68=="SOMETIMES")
	
strb_69_num <- rep(NA,n)
strb_69_num[which(cchs_mod$strb_69=="OFTEN")] <- 3
strb_69_num[which(cchs_mod$strb_69=="SOMETIMES")] <- 2
strb_69_num[which(cchs_mod$strb_69=="RARELY")] <- 1
strb_69_num[which(cchs_mod$strb_69=="NEVER")] <- 0

strb_69_OFTEN 	  <- (cchs_mod$strb_69=="OFTEN")
strb_69_SOMETIMES <- (cchs_mod$strb_69=="SOMETIMES")
strb_69_RARELY    <- (cchs_mod$strb_69=="RARELY")
strb_69_YES       <- (cchs_mod$strb_69=="OFTEN")|(cchs_mod$strb_69=="SOMETIMES")
	
dhhb_sex_num <- rep(NA,n)
dhhb_sex_num[which(cchs_mod$dhhb_sex=="MALE")] <- 1
dhhb_sex_num[which(cchs_mod$dhhb_sex=="FEMALE")] <- 0
	
ssmbdtng_num <- cchs_mod$ssmbdtng
ssmbdtng_num[ssmbdtng_num==99] <- NA

ssmbdtng_grt_q1 <- (cchs_mod$ssmbdtng>quantile(cchs_mod$ssmbdtng,0.25))
ssmbdtng_lst_q1 <- (cchs_mod$ssmbdtng<quantile(cchs_mod$ssmbdtng,0.25))
ssmbdtng_grt_q2 <- (cchs_mod$ssmbdtng>quantile(cchs_mod$ssmbdtng,0.5))
ssmbdtng_lst_q2 <- (cchs_mod$ssmbdtng<quantile(cchs_mod$ssmbdtng,0.5))
ssmbdtng_grt_q3 <- (cchs_mod$ssmbdtng>quantile(cchs_mod$ssmbdtng,0.75))
ssmbdtng_lst_q3 <- (cchs_mod$ssmbdtng<quantile(cchs_mod$ssmbdtng,0.75))
			
ssmbdsoc_num <- cchs_mod$ssmbdsoc
ssmbdsoc_num[ssmbdsoc_num==99] <- NA

ssmbdsoc_grt_q1 <- (cchs_mod$ssmbdsoc>quantile(cchs_mod$ssmbdsoc,0.25))
ssmbdsoc_lst_q1 <- (cchs_mod$ssmbdsoc<quantile(cchs_mod$ssmbdsoc,0.25))
ssmbdsoc_grt_q2 <- (cchs_mod$ssmbdsoc>quantile(cchs_mod$ssmbdsoc,0.5))
ssmbdsoc_lst_q2 <- (cchs_mod$ssmbdsoc<quantile(cchs_mod$ssmbdsoc,0.5))
ssmbdsoc_grt_q3 <- (cchs_mod$ssmbdsoc>quantile(cchs_mod$ssmbdsoc,0.75))
ssmbdsoc_lst_q3 <- (cchs_mod$ssmbdsoc<quantile(cchs_mod$ssmbdsoc,0.75))
			
ssmbdaff_num <- cchs_mod$ssmbdaff
ssmbdaff_num[ssmbdaff_num==99] <- NA

ssmbdaff_grt_q1 <- (cchs_mod$ssmbdaff>quantile(cchs_mod$ssmbdaff,0.25))
ssmbdaff_lst_q1 <- (cchs_mod$ssmbdaff<quantile(cchs_mod$ssmbdaff,0.25))
ssmbdaff_grt_q2 <- (cchs_mod$ssmbdaff>quantile(cchs_mod$ssmbdaff,0.5))
ssmbdaff_lst_q2 <- (cchs_mod$ssmbdaff<quantile(cchs_mod$ssmbdaff,0.5))
ssmbdaff_grt_q3 <- (cchs_mod$ssmbdaff>quantile(cchs_mod$ssmbdaff,0.75))
ssmbdaff_lst_q3 <- (cchs_mod$ssmbdaff<quantile(cchs_mod$ssmbdaff,0.75))
			
ssmbdemo_num <- cchs_mod$ssmbdemo
ssmbdemo_num[ssmbdemo_num==99] <- NA

ssmbdemo_grt_q1 <- (cchs_mod$ssmbdemo>quantile(cchs_mod$ssmbdemo,0.25))
ssmbdemo_lst_q1 <- (cchs_mod$ssmbdemo<quantile(cchs_mod$ssmbdemo,0.25))
ssmbdemo_grt_q2 <- (cchs_mod$ssmbdemo>quantile(cchs_mod$ssmbdemo,0.5))
ssmbdemo_lst_q2 <- (cchs_mod$ssmbdemo<quantile(cchs_mod$ssmbdemo,0.5))
ssmbdemo_grt_q3 <- (cchs_mod$ssmbdemo>quantile(cchs_mod$ssmbdemo,0.75))
ssmbdemo_lst_q3 <- (cchs_mod$ssmbdemo<quantile(cchs_mod$ssmbdemo,0.75))	
	
dhhbgage_num <- rep(NA,n)
dhhbgage_num[which(cchs_mod$dhhbgage=="15 TO 19 YEARS")] <- 17
dhhbgage_num[which(cchs_mod$dhhbgage=="20 TO 24 YEARS")] <- 22
dhhbgage_num[which(cchs_mod$dhhbgage=="25 TO 29 YEARS")] <- 27
dhhbgage_num[which(cchs_mod$dhhbgage=="30 TO 34 YEARS")] <- 32
dhhbgage_num[which(cchs_mod$dhhbgage=="35 TO 39 YEARS")] <- 37
dhhbgage_num[which(cchs_mod$dhhbgage=="40 TO 44 YEARS")] <- 42
dhhbgage_num[which(cchs_mod$dhhbgage=="45 TO 49 YEARS")] <- 47
dhhbgage_num[which(cchs_mod$dhhbgage=="50 TO 54 YEARS")] <- 52
dhhbgage_num[which(cchs_mod$dhhbgage=="55 TO 59 YEARS")] <- 57
dhhbgage_num[which(cchs_mod$dhhbgage=="60 TO 64 YEARS")] <- 62
dhhbgage_num[which(cchs_mod$dhhbgage=="65 TO 69 YEARS")] <- 67
dhhbgage_num[which(cchs_mod$dhhbgage=="70 TO 74 YEARS")] <- 72
dhhbgage_num[which(cchs_mod$dhhbgage=="75 TO 79 YEARS")] <- 77
dhhbgage_num[which(cchs_mod$dhhbgage=="80 YEARS OR MORE")] <- 82
	
dhhbgage1  <- cchs_mod$dhhbgage=="20 TO 24 YEARS"
dhhbgage2  <- cchs_mod$dhhbgage=="25 TO 29 YEARS"
dhhbgage3  <- cchs_mod$dhhbgage=="30 TO 34 YEARS"
dhhbgage4  <- cchs_mod$dhhbgage=="35 TO 39 YEARS"
dhhbgage5  <- cchs_mod$dhhbgage=="40 TO 44 YEARS"
dhhbgage6  <- cchs_mod$dhhbgage=="45 TO 49 YEARS"
dhhbgage7  <- cchs_mod$dhhbgage=="50 TO 54 YEARS"
dhhbgage8  <- cchs_mod$dhhbgage=="55 TO 59 YEARS"
dhhbgage9  <- cchs_mod$dhhbgage=="60 TO 64 YEARS"
dhhbgage10 <- cchs_mod$dhhbgage=="75 TO 79 YEARS"

dhhbgage_grt_20 <- (dhhbgage_num>20)
dhhbgage_lst_20 <- (dhhbgage_num<20)
dhhbgage_grt_25 <- (dhhbgage_num>25)
dhhbgage_lst_25 <- (dhhbgage_num<25)
dhhbgage_grt_30 <- (dhhbgage_num>30)
dhhbgage_lst_30 <- (dhhbgage_num<30)
dhhbgage_grt_35 <- (dhhbgage_num>35)
dhhbgage_lst_35 <- (dhhbgage_num<35)
dhhbgage_grt_40 <- (dhhbgage_num>40)
dhhbgage_lst_40 <- (dhhbgage_num<40)
			
edubdr04_1 <- cchs_mod$edubdr04=="SECONDARY GRAD."
edubdr04_2 <- cchs_mod$edubdr04=="OTHER POST-SEC."
edubdr04_3 <- cchs_mod$edubdr04=="POST-SEC. GRAD."	
edubdr04_4 <- (cchs_mod$edubdr04=="POST-SEC. GRAD.")|(cchs_mod$edubdr04=="OTHER POST-SEC.")
		
lbfbg31_num <- !((cchs_mod$lbfbg31=="NOT APPLICABLE")|(cchs_mod$lbfbg31=="NOT STATED"))
				
####################################################################################################

# Small design matrix so that we can see what 
# the "best" fit is by enumerating all subsets
# of models. Also for fitting trees
X_basic <- cbind(dhhb_sex_num,
				dhhbgage_num,
				edubdr04_1,					
				edubdr04_4,
				lbfbg31_num,
				strb_61_num,
				strb_62_num,
				strb_63_num,
				strb_64_num,
				strb_65A_num,
				strb_65B_num,
				strb_65C_num,
				strb_65D_num,
				strb_66_num,
				strb_67_num,
				strb_68_num,
				strb_69_num,
				ssmbdtng_num,
				ssmbdsoc_num,
				ssmbdaff_num,
				ssmbdemo_num)

colNames_basic <- c("sex",
					"age",
					"secGrad",
					"postGrad",
					"employed",
					"solve",
					"talk",
					"avoid",
					"sleep",
					"eat",
					"smoke",
					"drink",
					"drugMed",
					"exercise",
					"praySptHelp",
					"doEnjoyable",
					"brightSide",
					"tangSup",
					"socSup",
					"affSup",
					"emoSup")	
					
colnames(X_basic) <- colNames_basic


################################################################################

X_trans	<- cbind(
	dhhb_sex_num,
	dhhbgage_num,
	strb_61_num, 
	strb_61_OFTEN,
	strb_61_SOMETIMES,
	strb_61_RARELY,
	strb_61_YES,
	strb_62_num, 
	strb_62_OFTEN,
	strb_62_SOMETIMES,
	strb_62_RARELY,
	strb_62_YES,
	strb_63_num,
	strb_63_OFTEN,
	strb_63_SOMETIMES,
	strb_63_RARELY,
	strb_63_YES,
	strb_64_num,
	strb_64_OFTEN,
	strb_64_SOMETIMES,
	strb_64_RARELY,
	strb_64_YES,
	strb_65A_num,
	strb_65A_OFTEN,
	strb_65A_SOMETIMES,
	strb_65A_RARELY,
	strb_65A_YES,
	strb_65B_num,
	strb_65B_OFTEN,
	strb_65B_SOMETIMES,
	strb_65B_RARELY,
	strb_65B_YES,
	strb_65C_num,
	strb_65C_OFTEN,
	strb_65C_SOMETIMES,
	strb_65C_RARELY,
	strb_65C_YES,
	strb_65D_num,
	strb_65D_OFTEN,
	strb_65D_SOMETIMES,
	strb_65D_RARELY,
	strb_65D_YES,
	strb_66_num,
	strb_66_OFTEN,
	strb_66_SOMETIMES,
	strb_66_RARELY,
	strb_66_YES,
	strb_67_num,
	strb_67_OFTEN,
	strb_67_SOMETIMES,
	strb_67_RARELY,
	strb_67_YES,
	strb_68_num, 
	strb_68_OFTEN,
	strb_68_SOMETIMES,
	strb_68_RARELY,
	strb_68_YES,
	strb_69_num,
	strb_69_OFTEN,
	strb_69_SOMETIMES,
	strb_69_RARELY,
	strb_69_YES,
	ssmbdtng_num,
	ssmbdtng_grt_q1,
	ssmbdtng_lst_q1,
	ssmbdtng_grt_q2,
	ssmbdtng_lst_q2,
	ssmbdtng_grt_q3,
	ssmbdtng_lst_q3,		
	ssmbdsoc_num,
	ssmbdsoc_grt_q1,
	ssmbdsoc_lst_q1,
	ssmbdsoc_grt_q2,
	ssmbdsoc_lst_q2,
	ssmbdsoc_grt_q3,
	ssmbdsoc_lst_q3,	
	ssmbdaff_num,
	ssmbdaff_grt_q1,
	ssmbdaff_lst_q1,
	ssmbdaff_grt_q2,
	ssmbdaff_lst_q2,
	ssmbdaff_grt_q3,
	ssmbdaff_lst_q3,	
	ssmbdemo_num,
	ssmbdemo_grt_q1,
	ssmbdemo_lst_q1,
	ssmbdemo_grt_q2,
	ssmbdemo_lst_q2,
	ssmbdemo_grt_q3,
	ssmbdemo_lst_q3,
	dhhbgage1,
	dhhbgage2,
	dhhbgage3,
	dhhbgage4,
	dhhbgage5,
	dhhbgage6,
	dhhbgage_grt_20,
	dhhbgage_lst_20,
	dhhbgage_grt_25,
	dhhbgage_lst_25,
	dhhbgage_grt_30,
	dhhbgage_lst_30,
	dhhbgage_grt_35,
	dhhbgage_lst_35,
	dhhbgage_grt_40,
	dhhbgage_lst_40,
	edubdr04_1,
	edubdr04_2,
	edubdr04_3,
	edubdr04_4,
	lbfbg31_num)							

colNames_trans <- c(
	"sex",
	"age",
	"solve", 
	"solveOFTEN",
	"solveSOMETIMES",
	"solveRARELY",
	"solveYES",	
	"talk", 
	"talkOFTEN",
	"talkSOMETIMES",
	"talkRARELY",
	"talkYES",	
	"avoid", 
	"avoidOFTEN",
	"avoidSOMETIMES",
	"avoidRARELY",
	"avoidYES",	
	"sleep", 
	"sleepOFTEN",
	"sleepSOMETIMES",
	"sleepRARELY",
	"sleepYES",	
	"eat", 
	"eatOFTEN",
	"eatSOMETIMES",
	"eatRARELY",
	"eatYES",	
	"smoke", 
	"smokeOFTEN",
	"smokeSOMETIMES",
	"smokeRARELY",
	"smokeYES",	
	"drink", 
	"drinkOFTEN",
	"drinkSOMETIMES",
	"drinkRARELY",
	"drinkYES",		
	"drugMed", 
	"drugMedOFTEN",
	"drugMedSOMETIMES",
	"drugMedRARELY",
	"drugMedYES",	
	"exercise", 
	"exerciseOFTEN",
	"exerciseSOMETIMES",
	"exerciseRARELY",
	"exerciseYES",	
	"pray", 
	"prayOFTEN",
	"praySOMETIMES",
	"prayRARELY",
	"prayYES",		
	"enjoyable", 
	"enjoyableOFTEN",
	"enjoyableSOMETIMES",
	"enjoyableRARELY",
	"enjoyableYES",		
	"brightSide", 
	"brightSideOFTEN",
	"brightSideSOMETIMES",
	"brightSideRARELY",
	"brightSideYES",	
	"tangSup",
	"tangSupGrtQ1",	
	"tangSupLstQ1",		
	"tangSupGrtQ2",	
	"tangSupLstQ2",
	"tangSupGrtQ3",	
	"tangSupLstQ3",	
	"socSup",
	"socSupGrtQ1",	
	"socSupLstQ1",		
	"socSupGrtQ2",	
	"socSupLstQ2",
	"socSupGrtQ3",	
	"socSupLstQ3",
	"affSup",
	"affSupGrtQ1",	
	"affSupLstQ1",		
	"affSupGrtQ2",	
	"affSupLstQ2",
	"affSupGrtQ3",	
	"affSupLstQ3",
	"emoSup",
	"emoSupGrtQ1",	
	"emoSupLstQ1",		
	"emoSupGrtQ2",	
	"emoSupLstQ2",
	"emoSupGrtQ3",	
	"emoSupLstQ3",			
	"age20to24",
	"age25to29",
	"age30to34",
	"age35to39",
	"age40to44",
	"age45to49",		
	"ageGrt20",
	"ageLst20",
	"ageGrt25",
	"ageLst25",
	"ageGrt30",
	"ageLst30",		
	"ageGrt35",
	"ageLst35",	
	"ageGrt40",
	"ageLst40",		
	"secGrad",
	"othPostSec",
	"postSec",
	"postGrad",
	"employed")
	
colnames(X_trans) <- colNames_trans


####################################################################################################

colNames_interactions <- colNames_trans
X_interactions <- X_trans
for (i in 1:ncol(X_trans)) {
	if (i!=ncol(X_trans)) {
		for (j in (i+1):ncol(X_trans)) {
			X_interactions <- cbind(X_interactions,X_trans[,i]*X_trans[,j])
			colNames_interactions <- c(colNames_interactions,paste(colNames_trans[i],"*",colNames_trans[j],sep=""))
		}
	}
}
colnames(X_interactions) <- colNames_interactions	

####################################################################################################

inds <- which( apply(X_basic,2,sum.na.rm) >= MIN_OBS_IN_CATEGORY )
X_basic <- X_basic[,inds]


inds <- which( apply(X_trans,2,sum.na.rm) >= MIN_OBS_IN_CATEGORY )
X_trans <- X_trans[,inds]
		

inds <- which( apply(X_interactions,2,sum.na.rm) >= MIN_OBS_IN_CATEGORY )
X_interactions <- X_interactions[,inds]

####################################################################################################

if (TRUE) 
{	
	y <- cchs_mod$aldbfint
	X <- X_basic
	y1 <- rep(NA,n)
	y1[y=="YES"] <- 1
	y1[(y=="NO")|(y=="NOT APPLICABLE")] <- 0
	ind <- apply(cbind(y1,X),1,any.na)
	y1 <- y1[!ind]
	X1 <- X[!ind,]
	
	####################################################################################################
	
	#Xy <- as.data.frame(cbind(X1,y1))
	#bestAIC <- bestglm(Xy, IC="AIC",family=binomial)
	#bestBIC <- bestglm(Xy, IC="BIC",family=binomial)
	#bestEBIC <- bestglm(Xy, IC="BICg",family=binomial)
	#bestBICq <- bestglm(Xy, IC="BICq",family=binomial)

	####################################################################################################
		
	dat <- data.frame(y1,X1)
	
	print(summary( stepwise(y1,X1,gamma0=rep(0,ncol(X1)),type="AIC",binomial)$best.glm))
	
	####################################################################################################		
 
	library(rpart)
	
	control <- rpart.control(minsplit = 30, minbucket = 10, cp = 0.01,
		maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, xval = 10,
		surrogatestyle = 0, maxdepth = 30)
	res.rpart <- rpart(y1~.,data=dat,control=control)			
	
	library(rpart.plot)
	 
	#pdf("bpAlcInt.pdf",width=12,height=8)
	rpart.plot(res.rpart,type=4,extra=1,main="Alcohol Interference and Biploar",cex.main=2)
	#rpart.plot(res.rpart,type=4,extra=1,main="Alcohol Interference and Depression",cex.main=2)
	#rpart.plot(res.rpart,type=4,extra=1,main="Alcohol Interference and Mood",cex.main=2)
	#rpart.plot(res.rpart,type=4,extra=1,main="Alcohol Interference",cex.main=2)
	#dev.off()
	
	ans <- readline()
	
	####################################################################################################		
 
	y <- cchs_mod$aldbfint
	X <- X_trans
	y1 <- rep(NA,n)
	y1[y=="YES"] <- 1
	y1[(y=="NO")|(y=="NOT APPLICABLE")] <- 0
	ind <- apply(cbind(y1,X),1,any.na)
	y1 <- y1[!ind]
	X1 <- X[!ind,]
	
	res.stepwise <- stepwise(y1,X1,gamma0=rep(0,ncol(X1)),type="AIC",binomial)

	print(summary(res.stepwise$best.glm))
	
	gamma0 <- res.stepwise$best.gamma
	
	####################################################################################################		
 
	y <- cchs_mod$aldbfint
	X <- X_interactions
	y1 <- rep(NA,n)
	y1[y=="YES"] <- 1
	y1[(y=="NO")|(y=="NOT APPLICABLE")] <- 0
	ind <- apply(cbind(y1,X),1,any.na)
	y1 <- y1[!ind]
	X1 <- X[!ind,]

	
	print(summary( stepwise(y1,X1,gamma0=c(gamma0,rep(0,ncol(X1)-length(gamma0))),type="AIC",binomial)$best.glm))
}


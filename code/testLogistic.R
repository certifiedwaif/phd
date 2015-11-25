#!/usr/bin/env Rscript
# testLogistic.R
library(Rcpp)
setwd("~/Downloads/john/")
sourceCpp(file="fitLogistic.cpp")
load("CPM_IJC_AMELIA.RData")

X.clin <- Cleaned.data$data

X.des <- 	
	varNames <- c(  "Male",                 # Yes
	"Age",                   # Yes
	"Tum_NumNodesInv",       # Yes
	"TumMetSize",           # Yes
	"TumExtranodalY",       # Yes
	"TumCellType1",         # Yes
	"TumCellSize1",         # Yes
	"TumNecrosis",          # Yes
	"TumPigment",           # Yes
	"TumBRAFmutY",          # Yes
	"TumNRASmutY",          # Yes
	"PrimSiteSunExpC",        # Yes
	"PrimStage",         # Yes
	"PrimNaevusPresent",   # Yes
	"PrimBreslow",       # Yes
	"PrimMitos",         # Yes
	"PrimClark",          # Yes
	"PrimRegress",                   # Yes
	"PrimUlc",    # Yes
	"PrimHistolSSM",      # Yes
	"PrimHistolNM"       # Yes
)



P <- length(varNames)

X.des <- matrix(NA,nrow(X.clin),P)

 
X.des[,1] = X.clin$Person_Sex=="Male"
X.des[,2] = X.clin$Age_Analysis

X.des[,3] = as.numeric(as.character(X.clin$Tum_NumNodesInv))


X.des[,4] = as.numeric(as.character(X.clin$Tum_MetSize))
X.des[,5] = X.clin$Tum_Extranodal=="yes"
X.des[X.clin$Tum_CellType=="1",6] = 1
X.des[X.clin$Tum_CellType=="0",6] = 0
 
# Tum_CellSize
X.des[X.clin$Tum_CellSize=="0",7]  = 0
X.des[X.clin$Tum_CellSize=="1",7]  = 1
 

# Tum_Necrosis
X.des[,8] = X.clin$Tum_Necrosis
	

X.des[X.clin$Tum_Pigment=="1",9] = 1
X.des[X.clin$Tum_Pigment=="0",9] = 0

X.des[X.clin$Tum_BRAFmut=="Y",10] = 1
X.des[X.clin$Tum_BRAFmut=="N",10] = 0

X.des[X.clin$Tum_NRASmut=="Y",11] = 1
X.des[X.clin$Tum_NRASmut=="N",11] = 0

X.des[X.clin$Prim_Site_SunExp=="Chronic",12] = 1
X.des[X.clin$Prim_Site_SunExp=="Intermittent",12] = 0


X.des[X.clin$Prim_Stage=="Stage I",13] = 1
X.des[X.clin$Prim_Stage=="Stage II",13] = 2
X.des[X.clin$Prim_Stage=="Stage III",13] = 3

X.des[X.clin$Prim_Naevus=="Absent",14] = 0
X.des[X.clin$Prim_Naevus!="Absent",14] = 1
X.des[X.clin$Prim_Naevus=="NK",14]     = NA

X.des[,15] = X.clin$Prim_Breslow

X.des[,16] = X.clin$Prim_Mitos

X.des[X.clin$Prim_Clark==2,17] = 0 
X.des[X.clin$Prim_Clark==3,17] = 1
X.des[X.clin$Prim_Clark==4,17] = 2
X.des[X.clin$Prim_Clark==5,17] = 3


X.des[X.clin$Prim_Regress=="Absent",18] = 0
X.des[X.clin$Prim_Regress=="Early",18]  =  1
X.des[X.clin$Prim_Regress=="Late",18]   = 1	

X.des[X.clin$Prim_Ulc=="Yes",19] = 1
X.des[X.clin$Prim_Ulc=="No",19] = 0

X.des[X.clin$Prim_Histol!="Not Known",20] <- 0
X.des[X.clin$Prim_Histol=="SSM",20] = 1
X.des[X.clin$Prim_Histol=="SSM with NM",20] = 1

X.des[X.clin$Prim_Histol!="Not Known",21] <- 0
X.des[X.clin$Prim_Histol=="NM",21] = 1
X.des[X.clin$Prim_Histol=="NM (Minimal deviated melanoma)",21] = 1
X.des[X.clin$Prim_Histol=="NM (Spitz Naevus)",21] = 1
X.des[X.clin$Prim_Histol=="SSM with NM",21] = 1


X.des[14,3] <- NA  # Metsize is equal to tumour id

# Age 
# Tum_NumNodesInv 
# TumMetSize
# TumNecrosis
# PrimBreslow 
# PrimMitos

LogAge    <- log(1 + X.des[,2])
SqrtAge   <- sqrt(1 + X.des[,2])
SquareAge <- X.des[,2]^2
CubeAge   <- X.des[,2]^3



LogTum_NumNodesInv    <- log(1 + X.des[,3])
SqrtTum_NumNodesInv   <- sqrt(1 + X.des[,3])
SquareTum_NumNodesInv <- X.des[,3]^2
CubeTum_NumNodesInv   <- X.des[,3]^3



LogTumMetSize    <- log(1 + X.des[,4])
SqrtTumMetSize   <- sqrt(1 + X.des[,4])
SquareTumMetSize <- X.des[,4]^2
CubeTumMetSize   <- X.des[,4]^3



LogTumNecrosis    <- log(1 + X.des[,8])
SqrtTumNecrosis   <- sqrt(1 + X.des[,8])
SquareTumNecrosis <- X.des[,8]^2
CubeTumNecrosis   <- X.des[,8]^3



LogPrimBreslow    <- log(1 + X.des[,15])
SqrtPrimBreslow   <- sqrt(1 + X.des[,15])
SquarePrimBreslow <- X.des[,15]^2
CubePrimBreslow   <- X.des[,15]^3



LogPrimMitos    <- log(1 + X.des[,16])
SqrtPrimMitos   <- sqrt(1 + X.des[,16])
SquarePrimMitos <- X.des[,16]^2
CubePrimMitos   <- X.des[,16]^3

continuousVariables <- c(2,3,4,8,15,16) + 1

transAge <- cbind(LogAge,SqrtAge,SquareAge,CubeAge)
transTum_NumNodesInv <- cbind(LogTum_NumNodesInv,SqrtTum_NumNodesInv,SquareTum_NumNodesInv,CubeTum_NumNodesInv)
transTumMetSize <- cbind(LogTumMetSize,SqrtTumMetSize,SquareTumMetSize,CubeTumMetSize)
transTumNecrosis <- cbind(LogTumNecrosis,SqrtTumNecrosis,SquareTumNecrosis,CubeTumNecrosis)
transPrimBreslow <- cbind(LogPrimBreslow,SqrtPrimBreslow,SquarePrimBreslow,CubePrimBreslow)
transPrimMitos <- cbind(LogPrimMitos,SqrtPrimMitos,SquarePrimMitos,CubePrimMitos)


X.des <- cbind(X.des) #,transAge,transTum_NumNodesInv,transTumMetSize,transTumNecrosis,transPrimBreslow,transPrimMitos)




colnames(X.des) <- c(varNames)

if (FALSE) {
c(
"LogAge",
"SqrtAge",
"SquareAge",
"CubeAge",
"LogTum_NumNodesInv",
"SqrtTum_NumNodesInv",
"SquareTum_NumNodesInv",
"CubeTum_NumNodesInv",
"LogTumMetSize",
"SqrtTumMetSize",
"SquareTumMetSize",
"CubeTumMetSize",
"LogTumNecrosis",
"SqrtTumNecrosis",
"SquareTumNecrosis",
"CubeTumNecrosis",
"LogPrimBreslow",
"SqrtPrimBreslow",
"SquarePrimBreslow",
"CubePrimBreslow",
"LogPrimMitos",
"SqrtPrimMitos",
"SquarePrimMitos",
"CubePrimMitos")
}

 
vy <- as.numeric( Cleaned.data$data$Prognosis_status=="GP" )

 
mX <- cbind(1,X.des)
colnames(mX)[1] <- "Intercept"

n <- length(vy)
p <- ncol(mX)

source("fitLogistic.R")
 
#########################################################
		
# Imputation by the median value
for (j in 2:p) {
	vx <- mX[,j]
	if (any(is.na(mX[,j]))) {
		val <- median(vx,na.rm=TRUE)
		mX[is.na(vx),j] <- val
	}
}	

NUMVARS <- 6

models <- greycode(p-1)
vsize <- apply(models,1,sum) + 1
models <- models[vsize<=NUMVARS,]

print("Bohning")
library(lineprof)
a3 <- proc.time()[3];
res <- ALL_bohning(vy, mX, models, 1e-5, 1e-3) 
b3 <- proc.time()[3];   
cat(b3 - a3,"seconds \n") ;
cat("Models per second = ", dim(models)[1]/(b3 - a3), "\n");




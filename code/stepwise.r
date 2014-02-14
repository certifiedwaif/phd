
 
################################################################################

criteria <- function(y,X,gamma,type,family,maxterms)
{
    # Fit the model corresponding to gamma
    if (is.null(gamma)) {
    	res.glm <- glm(y~1,family=family) 
    } else {
    	if (sum(gamma)==0) {
    		res.glm <- glm(y~1,family=family) 
    	} else {
		    dat <- data.frame(y,X[,(gamma==1)])
		    res.glm <- glm(y~.,data=dat,family=family) 
    	}
	}
	
	if (any(is.na(res.glm$coef))) {
		score <- Inf
	} else {
		if (length(res.glm$coef)>=maxterms) {
			score <- Inf
		} else {
			if (any(summary(res.glm)$coef[-1,4]>0.5)) {
				score <- Inf
			} else {
				ll <- extractAIC(res.glm,k=0)
				if (type=="AIC") {
				    score <- ll[2] + 2*ll[1]	
				} 
				if (type=="BIC") {		
					score <- ll[2] + log(n)*ll[1]  
				}
			}
		}
	}   		
	return(list(score=score,res.glm=res.glm))
}

################################################################################

stepwise <- function(y,X,gamma0,type=c("AIC","BIC"),family,maxterms=sqrt(nrow(X))) 
{
    MAXMODELS <- 100000
    
	res <- criteria(y,X,gamma0,type,family,maxterms) 
	best.score <- res$score
	best.glm <- res$res.glm
	
	count <- 0
    
    FINISHED <- FALSE
    while(!FINISHED) 
    {
        CHANGED <- FALSE  
        
        # Flipping
        for (j in 1:ncol(X)) {
            count <- count + 1
			prop.gamma <- gamma0
			prop.gamma[j] <- 1 - prop.gamma[j]
			
			res <- criteria(y,X,prop.gamma,type,family,maxterms) 
			if (res$score < best.score) {
				cat("FLIPPED ",j,"\n")
				best.score <- res$score
				gamma0 <- prop.gamma
				best.glm <- res$res.glm
				CHANGED <- TRUE  
					
				cat("count=",count,"best.score=",best.score,"best.gamma=",which(gamma0==1),"\n")
				print(summary(res$res.glm))
				if (count > MAXMODELS) {
				   break;
				}
			}
        }
        
        # Swapping
        inds1 <- which(gamma0==1)    
        for (k1 in 1:length(inds1)) {
        	j1 <- inds1[k1]
	        inds0 <- which(gamma0==0)
        	for (k2 in 1:length(inds0)) {
        		j2 <- inds0[k2]
        		count <- count + 1
				prop.gamma <- gamma0
				prop.gamma[j1] <- 0
				prop.gamma[j2] <- 1
				res <- criteria(y,X,prop.gamma,type,family,maxterms) 
				if (res$score < best.score) {
					cat("SWAPPED ",j1," and ",j2,"\n")
					best.score <- res$score
					gamma0 <- prop.gamma
					best.glm <- res$res.glm
					CHANGED <- TRUE  	
					cat("count=",count,"best.score=",best.score,"best.gamma=",which(gamma0==1),"\n")
					print(summary(res$res.glm))
					if (count > MAXMODELS) {
					   break;
					}
				}
        	}
        	
        }
        
        if (!CHANGED) {
            break;
        }
    }

    return(list(best.gamma=gamma0,best.glm=best.glm))
}



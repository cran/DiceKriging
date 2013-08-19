`kmNoNugget.init` <-
function(model) {

	parinit <- model@parinit
	ninit <- model@control$pop.size
	param.n <- model@covariance@param.n      
	
	if (length(parinit)>0) {
	  matrixinit <- matrix(parinit, nrow = param.n, ncol = ninit) 
	} else {
	  lower <- model@lower
	  upper <- model@upper
	  if (existsMethod("paramSample", signature = class(model@covariance))) {
	    matrixinit <- paramSample(model@covariance, n=ninit, lower=lower, upper=upper, y=model@y)
	  } else {
	    # sample ninit design points, generated from uniform [lower, upper]
	    matrixinit <- matrix(runif(ninit*param.n), nrow = param.n, ncol = ninit)
	    matrixinit <- lower + matrixinit*(upper - lower)
	  }
		# take the best point				     
		logLikinit <- apply(matrixinit, 2, logLikFun, model)
		parinit <- matrixinit[, which.max(logLikinit)]    
	}
	
	model@parinit <- as.numeric(parinit)	
  
	return(model)
}

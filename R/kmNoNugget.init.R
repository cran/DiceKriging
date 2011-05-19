`kmNoNugget.init` <-
function(model) {

	parinit <- model@parinit
	
	if (identical(parinit, numeric(0))) {
		lower <- model@lower
		upper <- model@upper
		ninit <- model@control$pop.size
			# sample ninit design points, generated from uniform [lower, upper]
		param.n <- model@covariance@param.n
    matrixinit <- matrix(runif(ninit*param.n), param.n, ninit)
    
    if (!is(model@covariance, "covAffineScaling")) {
		  matrixinit <- lower + matrixinit*(upper - lower)
    } else {
      matrixinit <- 1/upper + matrixinit*(1/lower - 1/upper)
      matrixinit <- 1/matrixinit
    }
		
		# take the best point				     
		logLikinit <- apply(matrixinit, 2, logLikFun, model)
		parinit <- matrixinit[, which.max(logLikinit)]    
	}
	
	model@parinit <- as.numeric(parinit)	
	return(model)
}


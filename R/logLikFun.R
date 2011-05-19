`logLikFun` <-
function(param, model, envir=NULL) {
		
	if (identical(model@case, "NoNugget")) {
		
		model@covariance <- vect2covparam(param, model@covariance)
		model@covariance@sd2 <- 1		# to get the correlation matrix
		
		aux <- covMatrix(model@covariance, model@X)
		 
		R <- aux[[1]]
		T <- chol(R)
		    
   	x <- backsolve(t(T), model@y, upper.tri = FALSE)
   	M <- backsolve(t(T), model@F, upper.tri = FALSE)
		if (identical(model@known.param, "Trend")) {
   			z <- x - M %*% model@trend.coef
   		} else {
			Q <- qr.Q(qr(M))
   			H <- Q %*% t(Q)
   			z <- x - H %*% x
   		}
	
		sigma2.hat <- t(z)%*%(z) / model@n
		logLik <- -0.5*(model@n * log(2*pi*sigma2.hat) + 2*sum(log(diag(T))) + model@n)
		
		if (!is.null(envir)) { 
			assign("T",T, envir=envir) # stocking intermediate variables
			assign("R",R, envir=envir)
			assign("z",z, envir=envir)
			assign("sigma2.hat",sigma2.hat, envir=envir)
		}
		
	} 	else if (identical(model@case, "Nuggets")) {
		
		nparam <- length(param)
		
		model@covariance <- vect2covparam(param[1:(nparam-1)], model@covariance)
		model@covariance@sd2 <- param[nparam]
		
		aux <- covMatrix(model@covariance, model@X, noise.var=model@noise.var)
	
		C <- aux[[1]]
		C0 <- aux[[2]] 

    	T <- chol(C)
    	x <- backsolve(t(T), model@y, upper.tri = FALSE)
    	M <- backsolve(t(T), model@F, upper.tri = FALSE)
    	if (identical(model@known.param, "Trend")) {
   			z <- x - M %*% model@trend.coef
   		} else {
    		Q <- qr.Q(qr(M))
    		H <- Q %*% t(Q)
    		z <- x - H %*% x
		}
		logLik <-  -0.5*(model@n * log(2*pi) + 2*sum(log(diag(T))) + t(z)%*%z)     
	
		if (!is.null(envir)) {
			assign("T",T, envir=envir)	# stocking intermediate variables
			assign("C",C, envir=envir)
			assign("C0",C0, envir=envir)
			assign("z",z, envir=envir)
		}
				
	} else if (identical(model@case, "1Nugget")) {
	
		nparam <- length(param)
			
		model@covariance <- vect2covparam(param[1:(nparam-1)], model@covariance)
		model@covariance@sd2 <- 1
		model@covariance@nugget <- 0
		alpha <- param[nparam]
		
		aux <- covMatrix(model@covariance, model@X)
		R0 <- aux[[2]]
		R <- alpha*R0 + (1-alpha)*diag(model@n)
		#C0 <- aux[[2]] 

   		T <- chol(R)
   		x <- backsolve(t(T), model@y, upper.tri = FALSE)
   		M <- backsolve(t(T), model@F, upper.tri = FALSE)
   		if (identical(model@known.param, "Trend")) {
   			z <- x - M %*% model@trend.coef
   		} else {
	   		Q <- qr.Q(qr(M))
   			H <- Q %*% t(Q)
   			z <- x - H %*% x
		}
		v <- t(z)%*%(z) / model@n

		logLik <- -0.5*(model@n * log(2*pi*v) + 2*sum(log(diag(T))) + model@n)
	
		if (!is.null(envir)) {
			assign("T",T, envir=envir)	# stocking intermediate variables
			assign("R0",R0, envir=envir)
			assign("v",v, envir=envir)
			assign("z",z, envir=envir)
		}	
		
	}
	
	
	if (model@method=="PMLE") {
			fun <- match.fun(model@penalty$fun)
			penalty <- - model@n * sum(fun(1/model@covariance@range.val^2, model@penalty$value))
			logLik <- logLik + penalty
	}
		
	return(logLik)

}


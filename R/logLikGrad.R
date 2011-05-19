`logLikGrad` <-
function(param, model, envir) {
	
	if (identical(model@case, "NoNugget")) {
		
		toget <- matrix(c("R","T","z","sigma2.hat"),1,4)
		apply(toget, 2, get, envir=envir)

		model@covariance <- vect2covparam(param, model@covariance)
		model@covariance@sd2 <- 1		# to get the correlation matrix
	
		nparam <- length(param)
		
		logLik.derivative <- matrix(0,nparam,1)
									
		x <- backsolve(T,z)			# compute x := T^(-1)*z
		Rinv <- chol2inv(T)			# compute inv(R) by inverting T
			
		Rinv.upper <- Rinv[upper.tri(Rinv)]
		xx <- x%*%t(x)
		xx.upper <- xx[upper.tri(xx)]

		for (k in 1:nparam) {
			gradR.k <- covMatrixDerivative(model@covariance, X=model@X, C=R, k=k)
			gradR.k.upper <- gradR.k[upper.tri(gradR.k)]
		
			terme1 <- sum(xx.upper*gradR.k.upper)   / sigma2.hat   
				# quick computation of t(x)%*%gradR.k%*%x /  ...
			terme2 <- - sum(Rinv.upper*gradR.k.upper)       
				# quick computation of trace(Rinv%*%gradR.k)
			logLik.derivative[k] <- terme1 + terme2
		}
	
		if (model@method=="PMLE") {
			fun.derivative <- match.fun(model@penalty$fun.derivative)
			penalty <- - model@n * fun.derivative(1/param^2, model@penalty$value)*((-2)/param^3)
			logLik.derivative <- logLik.derivative + penalty
		}
	
	} else if (identical(model@case, "Nuggets")) {
	
		nparam <- length(param)
	
		model@covariance <- vect2covparam(param[1:(nparam-1)], model@covariance)
		model@covariance@sd2 <- sigma2 <- param[nparam]
		
		logLik.derivative <- matrix(0,nparam,1)
							
		toget <- matrix(c("C","T","C0","z"),1,4)
		apply(toget, 2, get, envir=envir)
	
		x <- backsolve(T,z)			# x := T^(-1)*z
		Cinv <- chol2inv(T)			# Invert R from given T
	
		Cinv.upper <- Cinv[upper.tri(Cinv)]
		xx <- x%*%t(x)
		xx.upper <- xx[upper.tri(xx)]

		# partial derivative with respect to parameters except sigma^2
		for (k in 1:(nparam-1)) {
			#gradC.k <- gradC[[k]]
			gradC.k <- covMatrixDerivative(model@covariance, X=model@X, C=C, k=k)
			gradC.k.upper <- gradC.k[upper.tri(gradC.k)]
			term1 <- sum(xx.upper*gradC.k.upper) #/ sigma2
				# economic computation of - t(x)%*%gradC.k%*%x / sigma2      
			term2 <- - sum(Cinv.upper*gradC.k.upper) 
				# economic computation of trace(Cinv%*%gradC.k)
			logLik.derivative[k] <- term1 + term2
		}
		# partial derivative with respect to sigma^2
		dCdsigma2 <- C0/sigma2
		term1 <- -t(x)%*%dCdsigma2%*%x 
		term2 <- sum(Cinv*dCdsigma2)       # economic computation of trace(Cinv%*%C0)
		logLik.derivative[nparam] <- -0.5*(term1 + term2) #/sigma2
	
		if (model@method=="PMLE") {
			fun.derivative <- match.fun(model@penalty$fun.derivative)
			param.pen <- model@covariance@range.val
			param.pen.pos <- 1:(nparam-1)
			penalty <- - model@n * fun.derivative(1/param.pen^2, model@penalty$value)*((-2)/param.pen^3)
			logLik.derivative[param.pen.pos] <- logLik.derivative[param.pen.pos] + penalty
		}
		
	} else if (identical(model@case, "1Nugget")) {
		
		nparam <- length(param)
	
		model@covariance <- vect2covparam(param[1:(nparam-1)], model@covariance)
		alpha <- param[nparam]
		model@covariance@sd2 <- 1  #sigma2 <- param[nparam-1]
		model@covariance@nugget <- 0  #param[nparam]
	
		logLik.derivative <- matrix(0,nparam,1)
							
		toget <- matrix(c("T","R0","v","z"),1,4)
		apply(toget, 2, get, envir=envir)
	
		x <- backsolve(T,z)			# x := T^(-1)*z
		Cinv <- chol2inv(T)			# Invert R from given T
	
		Cinv.upper <- Cinv[upper.tri(Cinv)]
		xx <- x%*%t(x)
		xx.upper <- xx[upper.tri(xx)]

		# partial derivative with respect to parameters except sigma^2
		for (k in 1:(nparam-1)) {
			gradC.k <- covMatrixDerivative(model@covariance, X=model@X, C=R0, k=k)
			gradC.k <- alpha*gradC.k
			gradC.k.upper <- gradC.k[upper.tri(gradC.k)]
			term1 <- sum(xx.upper*gradC.k.upper) / v
				# economic computation of - t(x)%*%gradC.k%*%x / v      
			term2 <- - sum(Cinv.upper*gradC.k.upper) 
				# economic computation of trace(Cinv%*%gradC.k)
			logLik.derivative[k] <- term1 + term2
		}
		# partial derivative with respect to v = sigma^2 + delta^2
		dCdv <- R0 - diag(model@n)
		term1 <- -t(x)%*%dCdv%*%x / v 
		term2 <- sum(Cinv*dCdv)       # economic computation of trace(Cinv%*%C0)
		logLik.derivative[nparam] <- -0.5*(term1 + term2) #/sigma2
	
		# partial derivative with respect to delta2
		#term1 <- -t(x)%*%x 
		#term2 <- sum(diag(Cinv))
		#logLik.derivative[nparam] <- -0.5*(term1 + term2)
	
		if (model@method=="PMLE") {
			fun.derivative <- match.fun(model@penalty$fun.derivative)
			param.pen <- model@covariance@range.val
			param.pen.pos <- 1:(nparam-2)
			penalty <- - model@n * fun.derivative(1/param.pen^2, model@penalty$value)*((-2)/param.pen^3)
			logLik.derivative[param.pen.pos] <- logLik.derivative[param.pen.pos] + penalty
		}

	}
	
	
	return(logLik.derivative)
}

